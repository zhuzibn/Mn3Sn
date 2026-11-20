function check_configuration(before_root)
%CHECK_CONFIGURATION CPU-only current interface checks; no external data needed.
% Supply before_root to additionally compare the historical capture.mat.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
history='not requested';
if nargin>0
    assert((ischar(before_root) && isrow(before_root) && ~isempty(before_root)) || ...
        (isstring(before_root) && isscalar(before_root) && ~ismissing(before_root) && strlength(before_root)>0), ...
        'Mn3Sn:ConfigurationFixture','before_root must be a nonempty historical capture directory.');
    file=fullfile(before_root,'capture.mat');
    assert(isfile(file),'Mn3Sn:ConfigurationFixture','Historical capture is missing: %s.',file);
    saved=load(file,'records','trajectories');
    assert(isfield(saved,'records') && iscell(saved.records) && numel(saved.records)==35 && ...
        isfield(saved,'trajectories') && iscell(saved.trajectories) && numel(saved.trajectories)==4, ...
        'Mn3Sn:ConfigurationFixture','Historical capture requires 35 preset records and four controls.');
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
    history='35 preset records and four open/short controls passed';
end
[presets,preset_failures]=check_presets;
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
fprintf('PASS configuration: %d presets, %d rejected preset calls, %d rejected configurations, effective edits and workspace isolation\n', ...
    presets,preset_failures,failures);
fprintf('Historical initialization comparison: %s\n',history);

    function reject(value)
        rejected=false;
        try,mn3sn_initialize(value);catch err,rejected=startsWith(err.identifier,'Mn3Sn:');end
        assert(rejected,'Invalid configuration accepted or failed outside validation.');
        failures=failures+1;
    end
end

function [count,failures]=check_presets
% Expected seed/current/field values exercise the public preset semantics.
cases={'relaxation',[],121,0,0;'fig2a',[],270,5e10,.01; ...
    'fig2b',[],270,5e10,-.01;'fig2c',[],270,-5e10,-.01;'fig4a',[],210,5e10,.01};
for name={'fig2d','fig5b','fig5c'}
    current=5e10;field=.01;
    if ~strcmp(name{1},'fig2d'),field=.1;end
    if strcmp(name{1},'fig5c'),current=-5e10;end
    for seed=30:60:330,cases(end+1,:)={name{1},seed,seed,current,field};end %#ok<AGROW>
end
for field=[175:25:400 450 500]
    cases(end+1,:)={'fig5a',field,210,5e10,field*1e-4}; %#ok<AGROW>
end
count=size(cases,1);failures=0;
for k=1:count
    row=cases(k,:);cfg=mn3sn_default_config(row{1},row{2});p=cfg.protocol;
    assert(strcmp(p.name,row{1}) && p.seed_deg==row{3} && p.current==row{4} && isequal(p.field,[0 row{5} 0]));
    assert(p.natomW==8 && p.natomL==8);
    if strcmp(p.name,'relaxation'),clock=[10e-12 .5e-15 1e-12 200 0 Inf];
    else,clock=[5e-9 5e-15 10e-12 50 1e-9 4e-9];end
    assert(isequal([p.runtime p.tstep p.gpusave p.savetstep p.pulse_start p.pulse_end],clock));
    assert(isequaln(cfg,mn3sn_validate_config(cfg)));
    state=mn3sn_initialize(cfg);assert(isequaln(state.config.protocol,p));
    assert(isequaln(p,mn3sn_protocol(row{1},row{2})));
    assert(isequaln(p,mn3sn_protocol(string(row{1}),row{2})));
    assert(isequaln(cfg,mn3sn_default_config(string(row{1}),row{2})));
end
assert(isequaln(mn3sn_default_config,mn3sn_default_config('relaxation')));
assert(isequaln(mn3sn_protocol,mn3sn_protocol('relaxation')));
for name={'relaxation','fig2a','fig2b','fig2c','fig2d','fig4a','fig5a','fig5b','fig5c'}
    value=[];
    if ismember(name{1},{'fig2d','fig5b','fig5c'}),value=210;
    elseif strcmp(name{1},'fig5a'),value=250;end
    cfg=mn3sn_default_config(name{1},value);
    assert(isequaln(cfg,mn3sn_default_config(name{1})) && isequaln(cfg.protocol,mn3sn_protocol(name{1})));
    for empty={[],single([]),uint16([]),sparse([]),''}
        assert(isequaln(cfg,mn3sn_default_config(name{1},empty{1})));
        assert(isequaln(cfg.protocol,mn3sn_protocol(name{1},empty{1})));
    end
end
for name={'fig2d','fig5b','fig5c','fig5a'}
    value=210;if strcmp(name{1},'fig5a'),value=300;end
    for invalid={single(value),uint16(value),sparse(value),complex(value,1),NaN,Inf, ...
            [value value],[value;value],true,'210',"210",{value},struct}
        reject({name{1},invalid{1}},'Mn3Sn:ProtocolValue');
    end
end
for name={'fig2d','fig5b','fig5c'},reject({name{1},31},'Mn3Sn:Seed');end
reject({'fig5a',280},'Mn3Sn:Field');
for name={'relaxation','fig2a','fig2b','fig2c','fig4a'},reject({name{1},30},'Mn3Sn:ProtocolValue');end
for name={'fig4b','RELAXATION',1,["fig2a" "fig2b"],['f';'i']},reject({name{1}},'Mn3Sn:Protocol');end

    function reject(args,identifier)
        for factory={@mn3sn_protocol,@mn3sn_default_config}
            call=factory{1};rejected=false;
            try,call(args{:});catch err,rejected=strcmp(err.identifier,identifier);end
            assert(rejected,'Preset input did not fail with %s.',identifier);
            failures=failures+1;
        end
    end
end
