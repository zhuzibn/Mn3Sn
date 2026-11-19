function result=benchmark_case(mode,name,root,execute,frozen)
%BENCHMARK_CASE Initialize or execute one production-path benchmark case.
if nargin<5,frozen=[];end
mode=validatestring(mode,{'quick','current'});
name=validatestring(name,{'relaxation','driven_nonuniform'});
selected_protocol=mn3sn_protocol('relaxation');
if strcmp(name,'driven_nonuniform')
    selected_protocol=mn3sn_protocol('fig4a');
    dt=selected_protocol.tstep;
    selected_protocol.runtime=60*dt;selected_protocol.gpusave=20*dt;
    selected_protocol.savetstep=1;
    selected_protocol.pulse_start=13*dt;selected_protocol.pulse_end=44*dt;
end
if strcmp(mode,'quick'),selected_protocol.natomW=4;selected_protocol.natomL=4;end
assert(selected_protocol.natomW==selected_protocol.natomL && ...
    selected_protocol.natomW==4+4*strcmp(mode,'current'),'Benchmark size changed.');
setup_timer=tic;
src=fileread(fullfile(root,'main.m'));
anchor='if(0)%view initial state';
assert(numel(strfind(src,anchor))==1,'Main initialization boundary changed.');
src=extractBefore(src,anchor);
% Validate dangerous options BEFORE evaluating any initialization code.
for token={'loadstartm=0;','load_fixed_atom_distrib=0;', ...
        'save_fixed_atom_distrib=0;','dwcalc=0;','thermalenable=0;', ...
        'dipolemode=0;','enablefixedge=0;','rk4=1;','bc=0;'}
    assert(numel(strfind(src,token{1}))==1,'Unsupported initialization: %s',token{1});
end
src=replace_once(src,'clear all;clc;close all;tic','');
src=replace_once(src,'gpuDevice(1)%select GPU device','');
lines=regexp(src,'(?m)^protocol=mn3sn_protocol[^\r\n]*','match');
assert(numel(lines)==1,'Protocol initialization changed.');
src=replace_once(src,lines{1},'protocol=selected_protocol;');
src=replace_once(src,"rng('shuffle');",'');
src=replace_once(src,"display('use random atom distribution')",'');
eval(src);
assert(rk4==1 && bc==0 && ~thermalenable && ~dipolemode && ~enablefixedge);
assert(~loadstartm && ~load_fixed_atom_distrib && ~save_fixed_atom_distrib && ~dwcalc);
assert(SOT_DLT==1 && SOT_FLT==0 && STT_DLT==0 && STT_FLT==0);
result=struct('schema_version',1,'mode',mode,'case_name',name);
result.atomtype=atomtype_;
result.axes=cat(3,coeffs_1,coeffs_2,coeffs_3);
result.generated_initial=cat(3,mx_init,my_init,mz_init);
% Snapshot effective production coefficients, controls and numerical settings.
names=strsplit(['gam ele mu_0 mub hbar kb rk4 bc DMIenable dwcalc thermalenable ' ...
    'loadstartm load_fixed_atom_distrib save_fixed_atom_distrib dipolemode enablefixedge ' ...
    'natomW natomL initial_theta_deg d natom Ksim Jintra Jfefe Jfegd gTM gamTM gamRE ' ...
    'tz volume musRE musTM msRE msTM Dsim alp jcSOT jcSTT Hext pulse_start pulse_end ' ...
    'SOT_DLT SOT_FLT psjSHE thetaSH chi BDSOTRE BDSOTTM STT_DLT STT_FLT psjSTT ' ...
    'etaSTT chiSTT BDSTTRE BDSTTTM T gpusave runtime tstep savetstep gpurun_number gpusteps totstep']);
result.config=struct;
for k=1:numel(names),result.config.(names{k})=eval(names{k});end
result.config.protocol=protocol;
occ=find(atomtype_==1);N=numel(occ);
assert(N==natomW*natomL*3/4 && all(ismember(atomtype_(:),[1 3])));
% Classes follow production initialization: a=seed, b=seed+120, c=seed+240.
[row,col]=ndgrid(1:natomW,1:natomL);classes=zeros(natomW,natomL);
classes(atomtype_==1 & mod(row,2)==1)=2;
classes(atomtype_==1 & ((mod(row,4)==2 & mod(col,2)==1) | ...
    (mod(row,4)==0 & mod(col,2)==0)))=3;
classes(atomtype_==1 & classes==0)=1;result.classes=classes;
assert(all(arrayfun(@(q)nnz(classes==q)==N/3,1:3)));
angles=(initial_theta_deg+(0:2)*120)*pi/180;
for q=1:3
    mask=classes==q;
    assert(max(abs(my_init(mask)-sin(angles(q))))<1e-14 && ...
        max(abs(mz_init(mask)-cos(angles(q))))<1e-14);
end
if ~isempty(frozen)
    for key={'schema_version','mode','case_name','atomtype','axes','classes','generated_initial','config'}
        assert(isequaln(result.(key{1}),frozen.(key{1})), ...
            'Frozen input configuration mismatch: %s',key{1});
    end
end
if strcmp(name,'driven_nonuniform')
    if isempty(frozen)
        % Algebraic fixture: deterministic, no RNG/version dependency.
        idx=(1:N)';m=[sin(idx),cos(sqrt(2)*idx),sin(sqrt(3)*idx+.3)];
        m=m./vecnorm(m,2,2);
        mx_init(occ)=m(:,1);my_init(occ)=m(:,2);mz_init(occ)=m(:,3);
    else
        mx_init=frozen.initial(:,:,1);my_init=frozen.initial(:,:,2);mz_init=frozen.initial(:,:,3);
    end
end
result.initial=cat(3,mx_init,my_init,mz_init);
if ~isempty(frozen),assert(isequal(result.initial,frozen.initial),'Frozen initial mismatch.');end
result.setup_seconds=toc(setup_timer);
if ~execute,return;end
dev=gpuDevice;wait(dev);solver_timer=tic;
integrate_llg;
wait(dev);result.integration_seconds=toc(solver_timer);
result.mmx=mmx;result.mmy=mmy;result.mmz=mmz;
result.saved_steps=saved_steps;result.t=t;
result.actual_updates=actual_updates;result.pulse_active_updates=pulse_active_updates;
result.observables=benchmark_validate(result);
end

function src=replace_once(src,old,new)
assert(numel(strfind(src,old))==1,'Initialization anchor changed: %s',old);
src=strrep(src,old,new);
end
