function check_coupled_rk4(out,phase)
% Execute the production loop in a local workspace; never run main.m/save.
% CPU checks replace only gpuArray transfers; GPU checks execute it unchanged.
if nargin<2, phase='initial'; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
S=load(fullfile(out,'static.mat'));
src=fileread(fullfile(out,'native','integrate_llg.m'));
assert(strcmp(src,fileread(fullfile(root,'integrate_llg.m'))),'Stale source snapshot');
src=strrep(src,'rknorm(atomtype_==3)=1;', ...
    'pre=max(pre,max(abs(rknorm(atomtype_==1)-1))); rknorm(atomtype_==3)=1;');
dev=gpuDevice(1);
if strcmp(phase,'initial')
    rng(20261102); torqueError=0;
    for q=1:40
        m=randn(1,3);m=m/norm(m);h=randn(1,3)*100;
        a=.003+.1*rand; ps=randn(1,3);ps=ps/norm(ps);pt=randn(1,3);pt=pt/norm(pt);
        b=randn(1,4);g=S.ref.gamma;sc=g/(1+a*a);
        [x,y,z]=atom_llgs_rhs(m(1),m(2),m(3),ps(1),ps(2),ps(3),pt(1),pt(2),pt(3), ...
            sc,a,h(1),h(2),h(3),b(1),b(2),b(3),b(4));r=[x y z];
        residual=r+g*cross(m,h)-a*cross(m,r) ...
            +g*b(1)*cross(m,cross(m,ps))-g*b(2)*cross(m,ps) ...
            +g*b(3)*cross(m,cross(m,pt))-g*b(4)*cross(m,pt);
        torqueError=max(torqueError,norm(residual)/norm(r));
    end
    assert(torqueError<1e-12);
    % Generic noncoplanar coupled trajectory: exposes frozen-field order loss.
    m=randn(3);m=m./vecnorm(m,2,2);p=S.native;
    dts=[10 5 2.5 1.25]*1e-15;states=cell(1,4);norms=zeros(1,4);
    for q=1:4
        [states{q},norms(q)]=simulate(src,p,m,round(.5e-12/dts(q)),dts(q),.01,p.b,false,2);
    end
    errors=zeros(1,3);
    for q=1:3,errors(q)=norm(states{q}(end,2:10)-states{q+1}(end,2:10));end
    order=log2(errors(1:2)./errors(2:3));assert(all(order>3.7 & order<4.3));
    [cpu,~]=simulate(src,p,m,20,5e-15,.01,p.b,false,2);
    [gpu,~]=simulate(src,p,m,20,5e-15,.01,p.b,true,2);
    gpuError=max(abs(cpu-gpu),[],'all');assert(gpuError<1e-11);
    % Larger source lattice, both legacy boundary flags: numerical invariants only.
    for bc=[0 1]
        [a,~]=simulate(src,S.native,[],8,.5e-15,0,0,true,5,bc);
        assert(all(isfinite(a),'all'));
    end
    writejson(fullfile(out,'production-gates.json'),struct('pass',true,'matlab',version, ...
        'gpu',dev.Name,'torque_relative_error',torqueError,'coupled_errors',errors, ...
        'coupled_order',order,'preprojection_errors',norms,'gpu_cpu_error',gpuError));
    fprintf('PASS production torque %.3g, coupled orders %.5f %.5f, GPU %.3g\n',torqueError,order,gpuError);
    cases=jsondecode(fileread(fullfile(out,'probe_cases.json')));
    cases=cases([cases.mode]==2);
else
    assert(strcmp(phase,'driven'));
    cases=jsondecode(fileread(fullfile(out,'production-replays.json')));
end
errors=zeros(1,numel(cases));
for q=1:numel(cases)
    c=cases(q);dt=c.dt_fs*1e-15;
    actual=readbin(fullfile(out,[c.name '.bin']));
    m=reshape(actual(1,2:10),3,3)';n=size(actual,1)-1;
    on=round(c.on_ns*1e-9/dt);off=on+round(c.width_ns*1e-9/dt);
    cuts=unique([0 on off n]);trajectory=actual(1,:);
    for j=1:numel(cuts)-1
        first=cuts(j);last=cuts(j+1);if last==first,continue;end
        drive=first>=on && first<off;
        [segment,~]=simulate(src,S.native,m,last-first,dt,c.field_Oe*1e-4*drive, ...
            S.native.b*c.current/5e10*drive,true,2);
        segment(:,1)=segment(:,1)+first*dt;
        trajectory=[trajectory;segment(2:end,:)]; %#ok<AGROW>
        m=reshape(segment(end,2:10),3,3)';
    end
    errors(q)=max(abs(trajectory(:,2:10)-actual(:,2:10)),[],'all');
    assert(max(abs(trajectory(:,1)-actual(:,1)))<1e-22);
    assert(errors(q)<1e-10,'Production GPU/replay mismatch: %s %.3g',c.name,errors(q));
    fprintf('PASS production GPU %s %.3g\n',c.name,errors(q));
end
writejson(fullfile(out,['production-' phase '.json']),struct('pass',true,'errors',errors,'cases',cases));
end

function [trajectory,pre]=simulate(src,p,m,n,dt,B,b,usegpu,W,bc)
if nargin<10,bc=1;end
assert(n>=3);constantfile;
natomW=W;natomL=W;dwcalc=0;load_fixed_atom_distrib=0;save_fixed_atom_distrib=0;
systemgeneration;
map=[4 1 2];
if W==2
    mx_init(map)=m(:,1);my_init(map)=m(:,2);mz_init(map)=m(:,3);
end
Jintra=p.J;Ksim=p.Kfield/2;Dsim=p.D;musRE=p.mu;gamRE=p.gamma;alp=p.alpha;
BDSOTRE=b;BDSTTRE=0;chi=0;dipolemode=0;thermalenable=0;enablefixedge=0;rk4=1;
psjSHEx=1;psjSHEy=0;psjSHEz=0;psjSTTx=0;psjSTTy=0;psjSTTz=1;
Hext=[0 B 0];tstep=dt;gpusteps=n+1;gpusave=n*dt;runtime=gpusave;
gpurun_number=1;savetstep=1;pre=0;
if ~usegpu
    src=strrep(src,'gpuArray(','double(');
    src=strrep(src,'arrayfun(','cpu_arrayfun(');
end
eval(src);
assert(size(mmx,3)==n+1);
assert(actual_updates==n && max(abs(t-(0:n)*dt))<1e-24);
mask=repmat(atomtype_==1,1,1,n+1);empty=repmat(atomtype_==3,1,1,n+1);
norms=sqrt(mmx.^2+mmy.^2+mmz.^2);assert(max(abs(norms(mask)-1))<1e-12);
assert(all(mmx(empty)==0 & mmy(empty)==0 & mmz(empty)==0));
assert(all(isfinite(mmx),'all') && all(isfinite(mmy),'all') && all(isfinite(mmz),'all'));
trajectory=zeros(n+1,10);trajectory(:,1)=(0:n)'*dt;
if W==2
    xx=reshape(mmx,4,[]);yy=reshape(mmy,4,[]);zz=reshape(mmz,4,[]);
    for j=1:3,trajectory(:,3*j-1:3*j+1)=[xx(map(j),:)' yy(map(j),:)' zz(map(j),:)'];end
end
end

function varargout=cpu_arrayfun(fun,varargin)
% GPU arrayfun expands scalar inputs; CPU arrayfun needs explicit expansion.
sz=size(varargin{1});
for k=1:numel(varargin)
    if isscalar(varargin{k}),varargin{k}=repmat(varargin{k},sz);end
end
[varargout{1:nargout}]=arrayfun(fun,varargin{:});
end

function data=readbin(path)
f=fopen(path,'rb','ieee-le');assert(f>=0);cleanup=onCleanup(@() fclose(f));data=fread(f,[10 Inf],'double')';
end
function writejson(path,data)
f=fopen(path,'w');assert(f>=0);cleanup=onCleanup(@() fclose(f));fprintf(f,'%s\n',jsonencode(data,PrettyPrint=true));
end
