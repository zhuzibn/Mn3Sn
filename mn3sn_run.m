function result=mn3sn_run(cfg,output_root)
%MN3SN_RUN Validated deterministic GPU LLGS with explicit, protected results.
% Omit output_root for a fresh external run directory; '' returns in memory.
root=fileparts(mfilename('fullpath'));
original=pwd;original_path=path;original_rng=rng;
cleanup=onCleanup(@() restore(original,original_path,original_rng));
% Bind all scripts/helpers to this source tree, independent of caller cwd/path.
cd(root);addpath(root,'-begin');
cfg=mn3sn_validate_config(cfg);
setup_timer=tic;state=mn3sn_initialize(cfg);
if nargin<2
    storage=fullfile(fileparts(root),'simulation_results','Mn3Sn');
    [~,unique_name]=fileparts(tempname);
    output_root=fullfile(storage,['run-' unique_name]);
elseif isequal(output_root,[])
    output_root='';
end
assert((ischar(output_root) && (isrow(output_root) || isempty(output_root))) || ...
    (isstring(output_root) && isscalar(output_root) && ~ismissing(output_root)), ...
    'Mn3Sn:OutputPath','Output must be a path, or empty for an in-memory run.');
output_root=char(output_root);
if ~isempty(output_root) && ~java.io.File(output_root).isAbsolute()
    output_root=fullfile(original,output_root);
end
if ~isempty(output_root),output_root=mn3sn_prepare_output(output_root);end
if isempty(cfg.rng_seed),rng('shuffle');else,rng(cfg.rng_seed,'twister');end
rng_before=rng;
% Indexed selection resets even an already-selected device. Preserve warm state
% and the caller's gpuArrays when the requested device is already current.
dev=gpuDevice;
if dev.Index~=cfg.gpu_device,dev=gpuDevice(cfg.gpu_device);end
provenance=source_provenance(root,dev);
setup_seconds=toc(setup_timer);
result=integrate(state,dev);
result.schema_name='Mn3Sn production result';result.schema_version=1;
result.requested_config=cfg;result.config=state.config;
result.atomtype=state.atomtype;result.axes=state.axes;result.classes=state.classes;
result.occupied_indices=find(state.atomtype==1);
result.generated_initial=state.generated_initial;result.initial=state.initial;
result.solver=struct('name','coupled_RK4','pulse_on_step',round(state.config.pulse_start/state.config.tstep), ...
    'pulse_off_step',round(state.config.pulse_end/state.config.tstep), ...
    'gate_convention','zero-based completed updates, on-inclusive/off-exclusive, constant across stages', ...
    'normalization','occupied spins after combined RK4 update only');
result.array_order='(W,L,time), MATLAB column-major site indexing';
result.site_meaning='1=occupied Mn, 3=vacant; classes a/b/c are initialization classes';
result.units=struct('spin','dimensionless','time','s','field_and_torque','T', ...
    'energy','J','moment','J/T','gamma','rad/(s T)','magnetization','A/m', ...
    'current_density','A/m^2','length','m','volume','m^3','seed_angle','deg', ...
    'temperature','K','hbar','eV s','sample_stride','updates');
result.output_root=output_root;result.setup_seconds=setup_seconds;
result.provenance=provenance;
result.provenance.rng_before=rng_before;result.provenance.rng_after=rng;
result.provenance.rng_seed=cfg.rng_seed;
result.provenance.rng_policy='[]: shuffle; integer: twister seed; caller RNG restored';
result.limitations={ ...
    'Class 3 code graph; physical crystallographic mapping and exchange-shell normalization unresolved', ...
    'Deterministic coupled RK4 only; no thermal, dipole, pinning, DW, restart, FL SOT or STT support', ...
    'Legacy open receiver weights are nonreciprocal; provenance does not establish physical validity'};
result.invariants=validate_result(result);
if ~isempty(output_root)
    file=fullfile(output_root,'result.mat');
    assert(~isfile(file),'Mn3Sn:OutputExists','Result file appeared during integration: %s.',file);
    save(file,'result','-v7.3');
end
end

function result=integrate(state,dev)
% Function-local aliases only; state fields are built from the validated whitelist.
names=fieldnames(state.config);
for k=1:numel(names),eval([names{k} '=state.config.(names{k});']);end
names=fieldnames(state.lattice);
for k=1:numel(names),eval([names{k} '=state.lattice.(names{k});']);end
atomtype_=state.atomtype;
coeffs_1=state.axes(:,:,1);coeffs_2=state.axes(:,:,2);coeffs_3=state.axes(:,:,3);
mx_init=state.initial(:,:,1);my_init=state.initial(:,:,2);mz_init=state.initial(:,:,3);
psjSHEx=psjSHE(1);psjSHEy=psjSHE(2);psjSHEz=psjSHE(3);
psjSTTx=psjSTT(1);psjSTTy=psjSTT(2);psjSTTz=psjSTT(3);
wait(dev);solver_timer=tic;
integrate_llg;
wait(dev);elapsed=toc(solver_timer);
result=struct('mmx',mmx,'mmy',mmy,'mmz',mmz,'t',t,'saved_steps',saved_steps, ...
    'actual_updates',actual_updates,'pulse_active_updates',pulse_active_updates, ...
    'integration_seconds',elapsed);
end

function invariants=validate_result(r)
c=r.config;steps=unique([0:c.savetstep:c.totstep c.totstep]);frames=numel(steps);
assert(r.actual_updates==c.totstep && isequal(r.saved_steps,steps) && ...
    isequal(size(r.t),size(steps)) && all(isfinite(r.t)) && ...
    max(abs(r.t-steps*c.tstep))<=8*eps(c.runtime),'Mn3Sn:Result','Result clock mismatch.');
on=round(c.pulse_start/c.tstep);off=round(c.pulse_end/c.tstep);
active=max(0,min(c.totstep,off)-min(c.totstep,on));
assert(r.pulse_active_updates==active,'Mn3Sn:Result','Pulse update count mismatch.');
occupied=repmat(r.atomtype==1,1,1,frames);empty=~occupied;
for a={r.mmx,r.mmy,r.mmz}
    assert(isa(a{1},'double') && isreal(a{1}) && all(isfinite(a{1}),'all') && ...
        isequal(size(a{1}),[c.natomW c.natomL frames]),'Mn3Sn:Result','Malformed trajectory.');
    assert(all(a{1}(empty)==0),'Mn3Sn:Result','Nonzero vacant spin.');
end
assert(isequal(cat(3,r.mmx(:,:,1),r.mmy(:,:,1),r.mmz(:,:,1)),r.initial), ...
    'Mn3Sn:Result','Initial frame mismatch.');
norms=sqrt(r.mmx.^2+r.mmy.^2+r.mmz.^2);error=max(abs(norms(occupied)-1));
assert(error<=1e-12,'Mn3Sn:Result','Occupied norm error exceeds 1e-12.');
invariants=struct('occupied_norm_error',error,'vacancy_max',0);
end

function p=source_provenance(root,dev)
p=struct('created_at',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss.SSSXXX')), ...
    'source_root',root,'matlab',version,'toolboxes',ver,'os',system_dependent('getos'), ...
    'gpu',struct('name',dev.Name,'index',dev.Index,'compute_capability',dev.ComputeCapability, ...
        'driver',dev.GraphicsDriverVersion,'toolkit',dev.ToolkitVersion,'total_memory',dev.TotalMemory));
p.git.root=git_query(root,'rev-parse --show-toplevel');
p.git.source_checkout=p.git.root.available && strcmpi( ...
    char(java.io.File(p.git.root.text).getCanonicalPath()),char(java.io.File(root).getCanonicalPath()));
p.git.reason='';
if ~p.git.source_checkout
    p.git.reason='Git is unavailable here, or its top-level directory is not this source tree.';
end
for row={{'head','rev-parse HEAD'},{'branch','branch --show-current'},{'status','status --short --branch'}}
    pair=row{1};
    if p.git.source_checkout,p.git.(pair{1})=git_query(root,pair{2});
    else,p.git.(pair{1})=struct('exit_status',[],'text','','available',false);end
end
files=dir(fullfile(root,'*.m'));
p.source=struct('name',{},'sha256',{},'bytes',{});
for k=1:numel(files)
    bytes=read_bytes(fullfile(root,files(k).name));
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
    hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    p.source(k)=struct('name',files(k).name,'sha256',hash,'bytes',bytes);
end
end

function value=git_query(root,args)
[status,output]=system(sprintf('git -C "%s" %s',root,args));
value=struct('exit_status',status,'text',strtrim(output),'available',status==0);
end

function bytes=read_bytes(file)
f=fopen(file,'rb');assert(f>=0,'Mn3Sn:Provenance','Cannot read source: %s.',file);
c=onCleanup(@() fclose(f));bytes=fread(f,Inf,'*uint8');
end

function restore(folder,oldpath,oldrng)
cd(folder);path(oldpath);rng(oldrng);
end
