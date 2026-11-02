function check_protocols(out)
% Actual main initialization and integration; no main output or result save.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);gpuDevice(1);
S=load(fullfile(out,'static.mat'));
presets={mn3sn_protocol('relaxation'),mn3sn_protocol('fig2a'), ...
    mn3sn_protocol('fig2b'),mn3sn_protocol('fig2c'),mn3sn_protocol('fig4a')};
for name={'fig2d','fig5b','fig5c'}
    for seed=30:60:330,presets{end+1}=mn3sn_protocol(name{1},seed);end %#ok<AGROW>
end
for field=[175:25:400 450 500],presets{end+1}=mn3sn_protocol('fig5a',field);end %#ok<AGROW>
for q=1:numel(presets)
    p=presets{q};r=source_run(p,false);
    assert(max(abs([r.mu/S.ref.mu,r.J/S.ref.J,r.Kfield/S.ref.Kfield, ...
        r.gamma/S.ref.gamma,r.Ms/S.ref.Ms,r.thickness/S.ref.thickness]-1))<1e-13);
    assert(abs(r.b-S.ref.b*p.current/5e10)<1e-17);
    assert(r.spacing==.4e-9 && r.thickness==30e-9);
    assert(r.gpusteps-1==round(p.gpusave/p.tstep));
    if strcmp(p.name,'relaxation')
        assert(r.occupied==16 && p.seed_deg==121 && p.current==0 && all(p.field==0));
        assert(r.updates==20000 && r.frames==101);
    else
        assert(r.occupied==3 && r.updates==1000000 && r.frames==20001);
        assert(p.pulse_start==1e-9 && p.pulse_end==4e-9);
        m=r.initial;angles=(p.seed_deg+(0:2)*120)'*pi/180;
        expected=[sin(angles)*cos(pi/2),sin(angles)*sin(pi/2),cos(angles)];
        assert(max(abs(m-expected),[],'all')<1e-14);
        switch p.name
            case 'fig2a',assert(p.current==5e10 && p.field(2)==.01 && p.seed_deg==270);
            case 'fig2b',assert(p.current==5e10 && p.field(2)==-.01 && p.seed_deg==270);
            case 'fig2c',assert(p.current==-5e10 && p.field(2)==-.01 && p.seed_deg==270);
            case 'fig2d',assert(p.current==5e10 && p.field(2)==.01);
            case 'fig4a',assert(p.current==5e10 && p.field(2)==.01 && p.seed_deg==210);
            case 'fig5a',assert(p.current==5e10 && p.seed_deg==210);
            case 'fig5b',assert(p.current==5e10 && p.field(2)==.1);
            case 'fig5c',assert(p.current==-5e10 && p.field(2)==.1);
        end
    end
end
for args={{'fig4b'},{'fig2d',31},{'fig5a',280},{'relaxation',30}}
    rejected=false;try,mn3sn_protocol(args{1}{:});catch,rejected=true;end
    assert(rejected,'Invalid protocol was accepted');
end
% Recompute an independent coupled triangle; pulses inside/on chunk edges.
errors=[];chunkErrors=[];
op=sum(S.operators(:,:,:,1),3);
for edges=[13 44;12 48]'
    for signs=[1 1;1 -1;-1 1;-1 -1]'
        p=mn3sn_protocol('fig4a');p.current=signs(1)*5e10;p.field(2)=signs(2)*.01;
        p.runtime=60*p.tstep;p.gpusave=p.runtime;p.savetstep=1;
        p.pulse_start=edges(1)*p.tstep;p.pulse_end=edges(2)*p.tstep;
        single=source_run(p,true);m=single.initial;expected=zeros(61,9);expected(1,:)=reshape(m',1,9);
        for k=0:59
            active=k>=edges(1)&&k<edges(2);
            field=@(u)reshape(op*u(:),3,3)+p.field*active;
            b=single.b*active;dt=p.tstep;
            k1=rhs(m,field(m),S.native,b);u=m+dt*k1/2;k2=rhs(u,field(u),S.native,b);
            u=m+dt*k2/2;k3=rhs(u,field(u),S.native,b);u=m+dt*k3;k4=rhs(u,field(u),S.native,b);
            m=m+dt/6*(k1+2*k2+2*k3+k4);m=m./vecnorm(m,2,2);
            expected(k+2,:)=reshape(m',1,9);
        end
        errors(end+1)=max(abs(single.spins-expected),[],'all'); %#ok<AGROW>
        p.gpusave=12*p.tstep;chunked=source_run(p,true);
        chunkErrors(end+1)=max(abs(single.spins-chunked.spins),[],'all'); %#ok<AGROW>
        assert(chunked.active_updates==diff(edges) && single.active_updates==diff(edges));
        assert(isequal(single.t,chunked.t) && numel(unique(chunked.t))==61);
        assert(max(abs(single.t-(0:60)*p.tstep))<1e-25);
        p.savetstep=7;sampled=source_run(p,true);indices=unique([0:7:60 60])+1;
        assert(isequal(sampled.spins,chunked.spins(indices,:)));
        assert(isequal(sampled.t,chunked.t(indices)));
    end
end
assert(max(errors)<1e-10 && max(chunkErrors)<1e-12);
% Actual default initialization, shortened duration only, on the 5x5 lattice.
p=mn3sn_protocol('relaxation');p.runtime=24*p.tstep;p.gpusave=8*p.tstep;p.savetstep=5;
source_run(p,true);
result=struct('pass',true,'presets_checked',numel(presets),'pulse_cases',numel(errors), ...
    'reference_errors',errors,'chunk_errors',chunkErrors,'matlab',version);
f=fopen(fullfile(out,'protocols.json'),'w');c=onCleanup(@() fclose(f));
fprintf(f,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('PASS %d presets, %d pulse cases; reference %.3g; chunk %.3g\n',numel(presets),numel(errors),max(errors),max(chunkErrors));
end

function r=source_run(selected_protocol,integrate)
root=fileparts(fileparts(mfilename('fullpath')));
src=fileread(fullfile(root,'main.m'));src=extractBefore(src,'if(0)%view initial state');
src=strrep(src,'clear all;clc;close all;tic','');src=regexprep(src,'gpuDevice\(1\)[^\n]*','');
src=regexprep(src,'protocol=mn3sn_protocol[^\n]*','protocol=selected_protocol;');
src=strrep(src,"display('use random atom distribution')",'');src=strrep(src,"rng('shuffle');",'');
eval(src);
r=struct('mu',musRE,'J',Jintra,'Kfield',2*Ksim,'gamma',gamRE,'Ms',msRE, ...
    'thickness',tz,'spacing',d,'b',BDSOTRE,'gpusteps',gpusteps,'updates',totstep, ...
    'occupied',nnz(atomtype_==1),'frames',numel(unique([0:savetstep:totstep totstep])));
map=[4 1 2];r.initial=[mx_init(map)' my_init(map)' mz_init(map)'];
if ~integrate,return;end
integrate_llg;
r.t=t;r.active_updates=pulse_active_updates;
assert(actual_updates==totstep && numel(t)==size(mmx,3));
assert(abs(t(end)-runtime)<1e-24 && t(1)==0 && all(diff(t)>0));
mask=repmat(atomtype_==1,1,1,numel(t));empty=repmat(atomtype_==3,1,1,numel(t));
norms=sqrt(mmx.^2+mmy.^2+mmz.^2);
assert(all(isfinite(norms),'all') && max(abs(norms(mask)-1))<1e-12);
assert(all(mmx(empty)==0 & mmy(empty)==0 & mmz(empty)==0));
if natomW==2
    x=reshape(mmx,4,[]);y=reshape(mmy,4,[]);z=reshape(mmz,4,[]);r.spins=zeros(numel(t),9);
    for i=1:3,r.spins(:,3*i-2:3*i)=[x(map(i),:)' y(map(i),:)' z(map(i),:)'];end
end
end

function r=rhs(m,h,p,b)
pol=repmat([1 0 0],3,1);
r=p.gamma/(1+p.alpha^2)*(-cross(m,h,2)-p.alpha*cross(m,cross(m,h,2),2) ...
    -b*cross(m,cross(m,pol,2),2)+p.alpha*b*cross(m,pol,2));
end
