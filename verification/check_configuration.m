function check_configuration(before_root)
%CHECK_CONFIGURATION CPU-only invalid-input and exact old/new initialization gates.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
saved=load(fullfile(before_root,'capture.mat'),'records','trajectories');
for k=1:numel(saved.records)
    old=saved.records{k};cfg=mn3sn_default_config;cfg.protocol=old.config.protocol;
    current=mn3sn_initialize(cfg);
    for key={'config','atomtype','axes','generated_initial'}
        assert(isequaln(current.(key{1}),old.(key{1})),'Initialization changed: %s, preset %d.',key{1},k);
    end
end
for k=1:numel(saved.trajectories)
    old=saved.trajectories{k};cfg=mn3sn_default_config;cfg.protocol=old.config.protocol;cfg.bc=old.config.bc;
    current=mn3sn_initialize(cfg);assert(isequaln(current.config,old.config));
    assert(isequal(current.atomtype,old.atomtype) && isequal(current.generated_initial,old.generated_initial));
end
cfg=mn3sn_default_config;failures=0;
bad=rmfield(cfg,'alp');reject(bad);bad=cfg;bad.alpp=.1;reject(bad);
bad=cfg;bad.protocol.extra=1;reject(bad);bad=cfg;bad.protocol=rmfield(bad.protocol,'tstep');reject(bad);
for value={NaN,Inf,1i,'0.1',[0 1],single(.1),sparse(.1)}
    bad=cfg;bad.alp=value{1};reject(bad);
end
for key={'rk4','thermalenable','dipolemode','enablefixedge','loadstartm', ...
        'load_fixed_atom_distrib','save_fixed_atom_distrib','dwcalc','SOT_FLT','STT_DLT','STT_FLT','chi','chiSTT','jcSTT'}
    bad=cfg;bad.(key{1})=1;if strcmp(key{1},'rk4'),bad.rk4=0;end;reject(bad);
end
bad=cfg;bad.rk4=2;reject(bad);
for row={{'natomW',5},{'natomL',3},{'natomW',4.5},{'savetstep',0},{'savetstep',1.5}, ...
        {'runtime',0},{'tstep',-1},{'gpusave',3e-16},{'runtime',10.1e-12}, ...
        {'pulse_start',-1},{'pulse_end',0},{'pulse_end',NaN},{'name','relax'}}
    v=row{1};bad=cfg;bad.protocol.(v{1})=v{2};reject(bad);
end
bad=cfg;bad.psjSHE=[0 0 0];reject(bad);bad=cfg;bad.psjSTT=[0;0;1];reject(bad);
bad=cfg;bad.protocol.field=[0 Inf 0];reject(bad);
for key={'d','gamRE','tz','volume','musRE','msRE','gpu_device'}
    bad=cfg;bad.(key{1})=0;reject(bad);
end
bad=cfg;bad.rng_seed=-1;reject(bad);bad=cfg;bad.rng_seed=2^32;reject(bad);
bad=cfg;bad.DMIenable=sparse(1);reject(bad);bad=cfg;bad.protocol.pulse_end=sparse(1);reject(bad);
bad=cfg;bad.msRE=realmin;bad.tz=realmin;reject(bad);
bad=cfg;bad.alp=realmax;reject(bad);bad=cfg;bad.Ksim=realmax;reject(bad);
bad=cfg;bad.initial_state=zeros(2,2,3);reject(bad);
initial=mn3sn_initialize(cfg);bad=cfg;bad.initial_state=initial.initial;
bad.initial_state(1,2,1)=1e-4;reject(bad);
bad=cfg;bad.initial_state=2*initial.initial;reject(bad);
bad=cfg;bad.initial_state=initial.initial;bad.initial_state(1)=NaN;reject(bad);
% Effective inputs must follow valid caller edits, not a hidden default block.
changed=cfg;changed.alp=.007;changed.musRE=cfg.musRE*1.1;
changed.msRE=cfg.msRE*1.2;changed.tz=cfg.tz*2;changed.protocol.current=-5e10;
changed.protocol.field=[.01 -.02 .03];changed.DMIenable=0;changed.rng_seed=7;
state=mn3sn_initialize(changed);c=state.config;
assert(c.alp==changed.alp && c.musRE==changed.musRE && c.Dsim==0 && isequal(c.Hext,changed.protocol.field));
assert(c.BDSOTRE==c.hbar/2*changed.thetaSH*changed.protocol.current/(changed.msRE*changed.tz));
changed.SOT_DLT=0;state=mn3sn_initialize(changed);assert(state.config.BDSOTRE==0);
changed=cfg;changed.protocol.pulse_start=1.4*cfg.protocol.tstep;
changed.protocol.pulse_end=5.6*cfg.protocol.tstep;mn3sn_validate_config(changed);
oldrng=rng;oldpath=path;oldpwd=pwd;
alp=99;mx_init=NaN;protocol=struct; %#ok<NASGU>
assert(isequaln(initial,mn3sn_initialize(cfg)) && isequal(oldrng,rng));
assert(strcmp(oldpwd,pwd) && strcmp(oldpath,path));
fprintf('PASS configuration: %d presets, %d open/short controls, %d rejected inputs, effective edits and workspace isolation\n', ...
    numel(saved.records),numel(saved.trajectories),failures);

    function reject(value)
        rejected=false;
        try,mn3sn_initialize(value);catch err,rejected=startsWith(err.identifier,'Mn3Sn:');end
        assert(rejected,'Invalid configuration accepted or failed outside validation.');
        failures=failures+1;
    end
end
