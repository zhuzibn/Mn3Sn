function check_three_spin_dynamics(out,stage)
% Source RHS, original CPU/GPU kernel, corrected reference, and replay gates.
S=load(fullfile(out,'static.mat')); assert(jsondecode(fileread(fullfile(out,'static.json'))).pass);
original=pwd; cleanup=onCleanup(@() cd(original));
addpath(fullfile(out,'native')); addpath(fullfile(out,'reference')); cd(out);
src=fileread(fullfile(out,'native','atomgpurk4.m'));
first=strfind(src,'u1=hhx;'); last=strfind(src,'kk1x=');
rhsExcerpt=src(first(1):last(1)-1);
rng(20261102); rhsError=0; gilbertError=0; oneStepError=0; gpuError=NaN;
for trial=1:20
    m=randn(3); m=m./vecnorm(m,2,2); h=randn(3)*100;
    for parameter=1:2
        if parameter==1, p=S.native; else, p=S.ref; end
        b=(-1)^trial*p.b; expected=rhs(m,h,p,b);
        for i=1:3
            got=source_rhs(m(i,:),h(i,:),p,b,rhsExcerpt);
            r=LLG_solver(p.alpha,m(i,:),h(i,:),[1 0 0],[0 0 0],[0 0 0],0,0,b,0)*p.gamma/(1+p.alpha^2);
            scale=max(1,norm(expected(i,:)));
            rhsError=max(rhsError,max([norm(got-expected(i,:)),norm(r-expected(i,:))])/scale);
            implicit=expected(i,:)+p.gamma*cross(m(i,:),h(i,:))-p.alpha*cross(m(i,:),expected(i,:)) ...
                +p.gamma*b*cross(m(i,:),cross(m(i,:),[1 0 0]));
            gilbertError=max(gilbertError,norm(implicit)/scale);
        end
        got=kernel(m,h,p,b,5e-15,false);
        expectedStep=diagnostic_step(m,p,b,5e-15,@(~)h);
        oneStepError=max(oneStepError,max(abs(got-expectedStep),[],'all'));
    end
end
assert(rhsError<1e-12&&gilbertError<1e-12&&oneStepError<1e-12);
% Verify native empty-slot reset on GPU with actual scalar kernel dispatch.
gpu=struct('available',false);
if gpuDeviceCount>0
    dev=gpuDevice(1); gpu=struct('available',true,'name',dev.Name,'compute_capability',dev.ComputeCapability);
    p=S.native; m=[S.geometry.initial([4 1 2],:);0 0 0];h=[ones(3)*.01;0 0 0];
    cpu=kernel(m,h,p,p.b,5e-15,false); got=kernel(m,h,p,p.b,5e-15,true);
    gpuError=max(abs(cpu-got),[],'all');assert(gpuError<1e-12&&all(got(4,:)==0));
end
fprintf('PASS RHS %.3g; Gilbert %.3g; frozen one-step %.3g; GPU %.3g\n',rhsError,gilbertError,oneStepError,gpuError);
% Short replay gate: all parameter/integration variants, current signs, dt and pulses.
if nargin<2,stage='initial';end
if strcmp(stage,'initial'),casefile='probe_cases.json';resultfile='dynamics.json';
else,assert(strcmp(stage,'driven'));casefile='driven_cases.json';resultfile='driven.json';end
cases=jsondecode(fileread(fullfile(out,casefile)));errors=zeros(numel(cases),1);
for c=1:numel(cases)
    spec=cases(c); if spec.mode==0, p=S.native; oi=1; else, p=S.ref; oi=2; end
    op=sum(S.operators(:,:,:,oi),3); dt=spec.dt_fs*1e-15;
    actual=readbin(fullfile(out,[spec.name '.bin']));
    m=reshape(actual(1,2:10),3,3)'; trajectory=zeros(size(actual));trajectory(1,:)=actual(1,:);
    n=round(spec.duration_ns*1e-9/dt); on=round(spec.on_ns*1e-9/dt);width=round(spec.width_ns*1e-9/dt);
    for k=0:n-1
        driven=k>=on&&k<on+width; B=[0 spec.field_Oe*1e-4*driven 0];b=p.b*spec.current/5e10*driven;
        field=@(u)reshape(op*u(:),3,3)+B;
        if spec.mode<2
            withEmpty=kernel([m;0 0 0],[field(m);0 0 0],p,b,dt,false);
            assert(all(withEmpty(4,:)==0));m=withEmpty(1:3,:);
        elseif spec.mode==2, m=diagnostic_step(m,p,b,dt,field);
        else, break; end
        trajectory(k+2,:)=[(k+1)*dt reshape(m',1,9)];
    end
    if spec.mode==3
        trajectory=reference_trajectory(out,spec,actual(1,2:10));
    end
    assert(isequal(size(trajectory),size(actual)));
    assert(max(abs(trajectory(:,1)-actual(:,1)))<1e-23);
    errors(c)=max(abs(trajectory(:,2:end)-actual(:,2:end)),[],'all');
    assert(errors(c)<1e-10,'Replay mismatch %s: %g',spec.name,errors(c));
    assert(all(isfinite(trajectory),'all'));
    fprintf('PASS %s replay max spin discrepancy %.3g\n',spec.name,errors(c));
end
writejson(fullfile(out,resultfile),struct('pass',true,'rhs_relative_error',rhsError, ...
    'gilbert_relative_error',gilbertError,'one_step_error',oneStepError,'gpu',gpu, ...
    'gpu_error',gpuError,'cases',cases,'replay_errors',errors));
end

function r=source_rhs(m,h,p,b,excerpt)
ssx=m(1);ssy=m(2);ssz=m(3);hhx=h(1);hhy=h(2);hhz=h(3);alph=p.alpha;
psjSHEx=1;psjSHEy=0;psjSHEz=0;psjSTTx=0;psjSTTy=0;psjSTTz=0;
BdSOT=b;BfSOT=0;BdSTT=0;BfSTT=0;eval(excerpt);
r=p.gamma/(1+p.alpha^2)*[dsdtx dsdty dsdtz];
end

function r=rhs(m,h,p,b)
pol=repmat([1 0 0],size(m,1),1);
r=p.gamma/(1+p.alpha^2)*(-cross(m,h,2)-p.alpha*cross(m,cross(m,h,2),2) ...
    -b*cross(m,cross(m,pol,2),2)+p.alpha*b*cross(m,pol,2));
end

function result=kernel(m,h,p,b,dt,usegpu)
if ~usegpu
    result=zeros(size(m));
    for i=1:size(m,1)
        [result(i,1),result(i,2),result(i,3)]=atomgpurk4(m(i,1),m(i,2),m(i,3),1,0,0,0,0,0, ...
            p.gamma/(1+p.alpha^2),p.alpha,dt,h(i,1),h(i,2),h(i,3),b,0,0,0);
    end
else
m=gpuArray(m);h=gpuArray(h);
[x,y,z]=arrayfun(@atomgpurk4,m(:,1),m(:,2),m(:,3),1,0,0,0,0,0, ...
    p.gamma/(1+p.alpha^2),p.alpha,dt,h(:,1),h(:,2),h(:,3),b,0,0,0);
result=gather([x y z]);
end
% Actual integrate_llg empty-slot dispatch policy, applied after normalization.
if size(m,1)==4,result(4,:)=0;end
end

function m=diagnostic_step(m,p,b,dt,field)
k1=rhs(m,field(m),p,b);u=m+dt*k1/2;k2=rhs(u,field(u),p,b);
u=m+dt*k2/2;k3=rhs(u,field(u),p,b);u=m+dt*k3;k4=rhs(u,field(u),p,b);
m=m+dt/6*(k1+2*k2+2*k3+k4);m=m./vecnorm(m,2,2);
end

function trajectory=reference_trajectory(out,spec,seed)
original=pwd; cleanup=onCleanup(@() cd(original));cd(fullfile(out,'reference'));initialize;
tstep=spec.dt_fs*1e-15;runtime=spec.duration_ns*1e-9;totstep=round(runtime/tstep)+1;
savestep=1;tstepcount=Inf;m=reshape(seed,3,3)';ma_init1=m(1,:);ma_init2=m(2,:);ma_init3=m(3,:);
on=round(spec.on_ns*1e-9/tstep);width=round(spec.width_ns*1e-9/tstep);
jc_SOTpulse=zeros(1,totstep);Hext_pulse=zeros(1,3,totstep);
jc_SOTpulse(on+1:on+width)=spec.current;
Hext_pulse(1,2,on+1:on+width)=spec.field_Oe*1e-4;
rk4_4llg;
trajectory=[(0:totstep-1)'*tstep,mmx(:,1),mmy(:,1),mmz(:,1), ...
    mmx(:,2),mmy(:,2),mmz(:,2),mmx(:,3),mmy(:,3),mmz(:,3)];
end

function data=readbin(path)
f=fopen(path,'rb','ieee-le');cleanup=onCleanup(@() fclose(f));data=fread(f,[10 Inf],'double')';
end
function writejson(path,data)
f=fopen(path,'w');cleanup=onCleanup(@() fclose(f));fprintf(f,'%s\n',jsonencode(data,PrettyPrint=true));
end
