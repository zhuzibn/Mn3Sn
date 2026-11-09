function check_periodic_production(out)
% Nonuniform production GPU trajectories against independently assembled fields.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);gpuDevice(1);
for dims=[5 5;4 2;6 8;8 7;0 8]'
    natomW=dims(1);natomL=dims(2);bc=0;
    rejected=false;
    try,distrib;catch err,rejected=strcmp(err.identifier,'Mn3Sn:PeriodicSize');end
    assert(rejected,'Invalid periodic dimensions were accepted.');
end
records={};
for W=[8 12 20]
    S=load(fullfile(out,sprintf('%dx%d-static.mat',W,W)));
    for signs=[1 1;1 -1;-1 1;-1 -1]'
        records{end+1}=replay(root,S,W,signs); %#ok<AGROW>
    end
end
result=struct('pass',true,'invalid_sizes_rejected',5,'cases',{[records{:}]},'matlab',version);
f=fopen(fullfile(out,'periodic-production.json'),'w');assert(f>=0);c=onCleanup(@() fclose(f));
fprintf(f,'%s\n',jsonencode(result,PrettyPrint=true));
fprintf('PASS periodic production: %d nonuniform driven GPU cases, maximum error %.3g\n', ...
    numel(records),max(cellfun(@(r)r.error,records)));
end

function result=replay(root,S,W,signs)
selected_protocol=mn3sn_protocol('fig4a');
selected_protocol.natomW=W;selected_protocol.natomL=W;
selected_protocol.current=signs(1)*5e10;selected_protocol.field=[0 signs(2)*.01 0];
dt=selected_protocol.tstep;
selected_protocol.runtime=60*dt;selected_protocol.gpusave=20*dt;selected_protocol.savetstep=1;
selected_protocol.pulse_start=13*dt;selected_protocol.pulse_end=44*dt;
src=fileread(fullfile(root,'main.m'));src=extractBefore(src,'if(0)%view initial state');
src=strrep(src,'clear all;clc;close all;tic','');src=regexprep(src,'gpuDevice\(1\)[^\n]*','');
src=regexprep(src,'protocol=mn3sn_protocol[^\n]*','protocol=selected_protocol;');
src=strrep(src,"display('use random atom distribution')",'');src=strrep(src,"rng('shuffle');",'');
eval(src);
occ=find(atomtype_==1);N=numel(occ);
assert(bc==0 && isequal(occ,S.occ) && isequal(atomtype_,S.atomtype_));
axes=[coeffs_1(occ) coeffs_2(occ) coeffs_3(occ)];assert(isequal(axes,S.axes));
rng(20261109);m=randn(N,3);m=m./vecnorm(m,2,2);
mx_init(occ)=m(:,1);my_init(occ)=m(:,2);mz_init(occ)=m(:,3);
expected=zeros(N,3,61);expected(:,:,1)=m;
for step=0:59
    active=step>=13 && step<44;
    field=@(u)-Jintra/musRE*S.A*u+Dsim/musRE*[zeros(N,1) S.C*u(:,3) -S.C*u(:,2)] ...
        +2*Ksim/musRE*sum(u.*axes,2).*axes+active*Hext;
    rhs=@(u)llg(u,field(u),gamRE,alp,BDSOTRE*active);
    k1=rhs(m);k2=rhs(m+dt*k1/2);k3=rhs(m+dt*k2/2);k4=rhs(m+dt*k3);
    m=m+dt/6*(k1+2*k2+2*k3+k4);m=m./vecnorm(m,2,2);
    expected(:,:,step+2)=m;
end
integrate_llg;
actual=zeros(N,3,61);
for k=1:61
    x=mmx(:,:,k);y=mmy(:,:,k);z=mmz(:,:,k);actual(:,:,k)=[x(occ) y(occ) z(occ)];
end
err=max(abs(actual-expected),[],'all');assert(err<1e-10);
assert(actual_updates==60 && pulse_active_updates==31 && isequal(saved_steps,0:60));
assert(max(abs(t-(0:60)*dt))<1e-25);
normError=max(abs(sqrt(sum(actual.^2,2))-1),[],'all');assert(normError<1e-12);
empty=repmat(atomtype_==3,1,1,61);
assert(all(mmx(empty)==0 & mmy(empty)==0 & mmz(empty)==0));
result=struct('W',W,'N',N,'current_sign',signs(1),'field_sign',signs(2), ...
    'error',err,'norm_error',normError,'updates',actual_updates,'pulse_updates',pulse_active_updates);
end

function r=llg(m,h,g,a,b)
p=repmat([1 0 0],size(m,1),1);
r=g/(1+a*a)*(-cross(m,h,2)-a*cross(m,cross(m,h,2),2) ...
    -b*cross(m,cross(m,p,2),2)+a*b*cross(m,p,2));
end
