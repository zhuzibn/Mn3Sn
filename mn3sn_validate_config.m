function cfg=mn3sn_validate_config(cfg)
%MN3SN_VALIDATE_CONFIG Validate the complete public input before any side effect.
% Physical/numerical inputs are real full doubles; flags also accept logicals.
root=fileparts(mfilename('fullpath'));original=pwd;original_path=path;
cleanup=onCleanup(@() restore(original,original_path));
cd(root);addpath(root,'-begin');
defaults=mn3sn_default_config;
fields_match(cfg,defaults,'cfg');
fields_match(cfg.protocol,defaults.protocol,'cfg.protocol');
p=cfg.protocol;
assert((ischar(p.name) && isrow(p.name)) || (isstring(p.name) && isscalar(p.name) && ~ismissing(p.name)), ...
    'Mn3Sn:Configuration','cfg.protocol.name must be a protocol name.');
p.name=char(p.name);
assert(ismember(p.name,{'relaxation','fig2a','fig2b','fig2c','fig2d','fig4a','fig5a','fig5b','fig5c'}), ...
    'Mn3Sn:Configuration','cfg.protocol.name is unknown.');
for key={'natomW','natomL','seed_deg','current','pulse_start','runtime','tstep','gpusave','savetstep'}
    scalar(p.(key{1}),['cfg.protocol.' key{1}]);
end
for key={'natomW','natomL'}
    v=p.(key{1});assert(v>=2 && v==fix(v) && v<=flintmax, ...
        'Mn3Sn:Configuration','cfg.protocol.%s must be an integer >=2.',key{1});
end
vector(p.field,'cfg.protocol.field',false);
assert(isa(p.pulse_end,'double') && isreal(p.pulse_end) && ~issparse(p.pulse_end) && isscalar(p.pulse_end) && ...
    ~isnan(p.pulse_end) && p.pulse_end>p.pulse_start && p.pulse_start>=0, ...
    'Mn3Sn:Configuration','Pulse edges require 0<=pulse_start<pulse_end (Inf end allowed).');
assert(p.runtime>0 && p.tstep>0 && p.gpusave>0 && p.savetstep>=1 && ...
    p.savetstep==fix(p.savetstep) && p.savetstep<=flintmax, ...
    'Mn3Sn:Configuration','Positive times and an integer savetstep>=1 are required.');
chunks=round(p.runtime/p.gpusave);steps=round(p.gpusave/p.tstep);
assert(isfinite(chunks*steps) && chunks>=1 && steps>=2 && chunks*steps<=flintmax && ...
    abs(chunks*p.gpusave-p.runtime)<1e-10*p.runtime && ...
    abs(steps*p.tstep-p.gpusave)<1e-10*p.gpusave && ...
    abs(chunks*steps*p.tstep-p.runtime)<1e-10*p.runtime, ...
    'Mn3Sn:TimeGrid','runtime/gpusave and gpusave/tstep must be integer counts, with >=2 updates per chunk.');
assert(p.pulse_start/p.tstep<=flintmax && ...
    (isinf(p.pulse_end) || p.pulse_end/p.tstep<=flintmax), ...
    'Mn3Sn:TimeGrid','Pulse indices must fit exact integer representation.');
flags={'bc','DMIenable','dwcalc','thermalenable','loadstartm','load_fixed_atom_distrib', ...
    'save_fixed_atom_distrib','enablefixedge','SOT_DLT','SOT_FLT','STT_DLT','STT_FLT'};
for key=flags
    v=cfg.(key{1});assert((isa(v,'double') || islogical(v)) && isreal(v) && ~issparse(v) && ...
        isscalar(v) && ismember(v,[0 1]),'Mn3Sn:Configuration','cfg.%s must be 0 or 1.',key{1});
    cfg.(key{1})=double(v);
end
for key={'rk4','dipolemode','d','Ksim','Jintra','gamRE','tz','volume','musRE','msRE', ...
        'Dsim','alp','thetaSH','chi','jcSTT','etaSTT','chiSTT','T','gpu_device'}
    scalar(cfg.(key{1}),['cfg.' key{1}]);
end
assert(cfg.rk4==1,'Mn3Sn:UnsupportedSolver','cfg.rk4 must be 1 (coupled RK4); legacy solvers are unsupported.');
for key={'dwcalc','thermalenable','loadstartm','load_fixed_atom_distrib','save_fixed_atom_distrib', ...
        'dipolemode','enablefixedge','SOT_FLT','STT_DLT','STT_FLT','chi','chiSTT','jcSTT'}
    assert(cfg.(key{1})==0,'Mn3Sn:UnsupportedFeature','cfg.%s must be 0 in the supported deterministic scope.',key{1});
end
if cfg.bc==0
    assert(p.natomW>=4 && mod(p.natomW,4)==0 && p.natomL>=4 && mod(p.natomL,2)==0, ...
        'Mn3Sn:PeriodicSize','Periodic arrays require W>=4 divisible by 4 and even L>=4.');
end
for key={'d','gamRE','tz','volume','musRE','msRE'}
    assert(cfg.(key{1})>0,'Mn3Sn:Configuration','cfg.%s must be positive.',key{1});
end
assert(cfg.alp>=0 && cfg.T>=0,'Mn3Sn:Configuration','cfg.alp and cfg.T must be nonnegative.');
assert(cfg.gpu_device>=1 && cfg.gpu_device==fix(cfg.gpu_device), ...
    'Mn3Sn:Configuration','cfg.gpu_device must be a positive integer.');
vector(cfg.psjSHE,'cfg.psjSHE',true);vector(cfg.psjSTT,'cfg.psjSTT',true);
if ~isequal(cfg.rng_seed,[])
    scalar(cfg.rng_seed,'cfg.rng_seed');
    assert(cfg.rng_seed>=0 && cfg.rng_seed<=2^32-1 && cfg.rng_seed==fix(cfg.rng_seed), ...
        'Mn3Sn:Configuration','cfg.rng_seed must be [] (shuffle) or an integer in [0,2^32-1].');
end
if ~isequal(cfg.initial_state,[])
    x=cfg.initial_state;
    assert(isa(x,'double') && isreal(x) && ~issparse(x) && ...
        isequal(size(x),[p.natomW p.natomL 3]) && all(isfinite(x),'all'), ...
        'Mn3Sn:InitialState','cfg.initial_state must be [] or finite real double (W,L,3).');
end
cfg.protocol=p;
end

function fields_match(value,template,label)
assert(isstruct(value) && isscalar(value),'Mn3Sn:Configuration','%s must be a scalar struct.',label);
missing=setdiff(fieldnames(template),fieldnames(value));extra=setdiff(fieldnames(value),fieldnames(template));
assert(isempty(missing),'Mn3Sn:Configuration','Missing %s fields: %s.',label,strjoin(missing,', '));
assert(isempty(extra),'Mn3Sn:Configuration','Unknown %s fields: %s.',label,strjoin(extra,', '));
end

function scalar(value,label)
assert(isa(value,'double') && isreal(value) && ~issparse(value) && isscalar(value) && isfinite(value), ...
    'Mn3Sn:Configuration','%s must be a finite real double scalar.',label);
end

function vector(value,label,unit)
assert(isa(value,'double') && isreal(value) && ~issparse(value) && ...
    isequal(size(value),[1 3]) && all(isfinite(value)), ...
    'Mn3Sn:Configuration','%s must be a finite real double 1x3 vector.',label);
if unit,assert(abs(norm(value)-1)<=1e-12,'Mn3Sn:Configuration','%s must have unit norm.',label);end
end

function restore(folder,oldpath)
cd(folder);path(oldpath);
end
