function check_production_interface(out,before_root)
%CHECK_PRODUCTION_INTERFACE Short GPU migration, schema and non-destructive I/O.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
out=mn3sn_prepare_output(out);
saved=load(fullfile(before_root,'capture.mat'),'trajectories');
original=pwd;original_path=path;original_rng=rng;
cleanup=onCleanup(@() restore(original,original_path,original_rng));
errors=zeros(1,numel(saved.trajectories));
for k=1:numel(saved.trajectories)
    old=saved.trajectories{k};cfg=mn3sn_default_config;
    cfg.protocol=old.config.protocol;cfg.bc=old.config.bc;cfg.rng_seed=7;
    result=mn3sn_run(cfg,'');
    for key={'config','atomtype','axes','generated_initial','saved_steps','actual_updates','pulse_active_updates'}
        assert(isequaln(result.(key{1}),old.(key{1})),'Migration changed %s.',key{1});
    end
    errors(k)=max([max(abs(result.mmx-old.mmx),[],'all'),max(abs(result.mmy-old.mmy),[],'all'), ...
        max(abs(result.mmz-old.mmz),[],'all')]);
    assert(errors(k)<=1e-10 && isequal(result.t,old.t));
end
assert(isequal(original_rng,rng) && strcmp(original,pwd) && strcmp(original_path,path));
% Persist a short case, then compare the saved schema byte-for-value to the return.
device_probe=gpuArray([1 2 3]);
destination=fullfile(out,'saved');result=mn3sn_run(cfg,destination);
assert(isequal(gather(device_probe),[1 2 3]),'Runner reset the already-selected GPU.');
loaded=load(fullfile(destination,'result.mat'));assert(isequal(fieldnames(loaded),{'result'}));
assert(isequaln(result,loaded.result));
assert(result.schema_version==1 && strcmp(result.schema_name,'Mn3Sn production result'));
assert(result.solver.pulse_on_step==round(cfg.protocol.pulse_start/cfg.protocol.tstep) && ...
    result.solver.pulse_off_step==round(cfg.protocol.pulse_end/cfg.protocol.tstep));
assert(isequal(result.requested_config,cfg) && result.provenance.rng_before.Seed==7);
assert(isequal(result.provenance.rng_before,result.provenance.rng_after));
assert(result.provenance.git.head.available && result.provenance.git.status.available);
assert(~isempty(result.provenance.created_at) && ~isempty(result.units.hbar));
for k=1:numel(result.provenance.source)
    item=result.provenance.source(k);
    assert(strcmp(item.sha256,hash_bytes(item.bytes)));
    assert(isequal(item.bytes,read_bytes(fullfile(root,item.name))));
end
oldmain=read_bytes(fullfile(before_root,'source','main.m'));
current_main=result.provenance.source(strcmp({result.provenance.source.name},'main.m'));
assert(~strcmp(hash_bytes(oldmain),current_main.sha256));
% No commit was needed to distinguish the old/new source content.
[status,head]=system(sprintf('git -C "%s" rev-parse HEAD',root));
assert(status==0 && strcmp(strtrim(head),result.provenance.git.head.text));
assert(isequal(original_rng,rng));
% Supported switches must reach the unchanged fields/torques during integration.
driven=mn3sn_default_config;driven.protocol=saved.trajectories{3}.config.protocol;
off=driven;off.DMIenable=0;zero=driven;zero.Dsim=0;
assert_same_trajectory(mn3sn_run(off,''),mn3sn_run(zero,''));
off=driven;off.SOT_DLT=0;zero=driven;zero.protocol.current=0;
assert_same_trajectory(mn3sn_run(off,''),mn3sn_run(zero,''));
% Exercise the actual thin wrapper with only its preset/output call substituted.
wrapped=run_wrapper(root,cfg);
assert(isequal(wrapped.mmx,result.mmx) && isequal(wrapped.mmy,result.mmy) && isequal(wrapped.mmz,result.mmz));
% Ordinary in-memory execution from an unrelated caller directory leaves it empty.
caller=fullfile(out,'caller');mkdir(caller);cd(caller);
for name={'constantfile.m','mn3sn_protocol.m'}
    f=fopen(fullfile(caller,name{1}),'w');assert(f>=0);
    fprintf(f,"error('Test:Shadow','Caller source must not supply production defaults');\n");fclose(f);
end
assert(isequaln(mn3sn_default_config,default_from_root(root)));
assert(isequaln(mn3sn_validate_config(cfg),cfg));
before_names={dir(caller).name};mn3sn_run(cfg,'');after_names={dir(caller).name};
assert(isequal(before_names,after_names) && strcmp(pwd,caller));
relative=mn3sn_run(cfg,'relative-run');
assert(strcmp(relative.output_root,fullfile(caller,'relative-run')));cd(original);

manifest_file=fullfile(root,'benchmark','baseline_manifest.json');
manifest=jsondecode(fileread(manifest_file));protected_files={manifest_file,fullfile(destination,'result.mat')};
for mode={'quick','current'}
    base=fullfile(root,manifest.(mode{1}).directory);
    for name={'relaxation.mat','driven_nonuniform.mat','metadata.json'}
        protected_files{end+1}=fullfile(base,name{1}); %#ok<AGROW>
    end
end
hashes=cellfun(@(p)hash_bytes(read_bytes(p)),protected_files,'UniformOutput',false);
rejections=0;baddevice=cfg;baddevice.gpu_device=1e6;
paths={destination,fullfile(destination,'result.mat'),fullfile(root,'forbidden-cq-output'), ...
    fullfile(out,'OneDrive','forbidden'),fullfile(root,manifest.quick.directory), ...
    fullfile(root,manifest.current.directory,'new-child'), ...
    fullfile(root,manifest.quick.directory,'unused','..','alias-child')};
for k=1:numel(paths)
    rejected=false;
    try,mn3sn_run(baddevice,paths{k});catch err,rejected=startsWith(err.identifier,'Mn3Sn:Output') || strcmp(err.identifier,'Mn3Sn:ProtectedBaseline');end
    assert(rejected,'Destination did not fail before selecting an invalid GPU: %s',paths{k});rejections=rejections+1;
end
for k=1:numel(protected_files),assert(strcmp(hashes{k},hash_bytes(read_bytes(protected_files{k}))));end
assert(~isfolder(fullfile(root,'forbidden-cq-output')) && ~isfolder(fullfile(out,'OneDrive')));
explicit=fullfile(out,'explicit-input');mkdir(explicit);
rejected=false;try,mn3sn_prepare_output(fullfile(explicit,'candidate'),explicit);catch err,rejected=strcmp(err.identifier,'Mn3Sn:ProtectedBaseline');end
assert(rejected && ~isfolder(fullfile(explicit,'candidate')));
claimed=mn3sn_prepare_output(fullfile(out,'competing-claim'));
rejected=false;try,mn3sn_prepare_output(claimed);catch err,rejected=strcmp(err.identifier,'Mn3Sn:OutputExists');end
assert(rejected);
% Initial-state validation happens before output directory creation.
bad=cfg;bad.initial_state=zeros(cfg.protocol.natomW,cfg.protocol.natomL,3);
invalid_out=fullfile(out,'invalid-initial');rejected=false;
try,mn3sn_run(bad,invalid_out);catch err,rejected=strcmp(err.identifier,'Mn3Sn:InitialState');end
assert(rejected && ~isfolder(invalid_out) && isequal(original_rng,rng));
% A source copy without Git or manifest must label metadata honestly and refuse persistence.
copy=fullfile(out,'source-without-git');mkdir(copy);
files=dir(fullfile(root,'*.m'));for k=1:numel(files),copyfile(fullfile(root,files(k).name),copy);end
cd(copy);addpath(copy,'-begin');
memory=mn3sn_run(cfg,'');
assert(~memory.provenance.git.head.available && ~memory.provenance.git.source_checkout && ...
    ~isempty(memory.provenance.git.reason));
rejected=false;missing_out=fullfile(out,'missing-manifest');
try,mn3sn_run(cfg,missing_out);catch err,rejected=strcmp(err.identifier,'Mn3Sn:MissingBaseline');end
assert(rejected && ~isfolder(missing_out));
cd(original);path(original_path);
summary=struct('passed',true,'trajectory_errors',errors,'save_load',true, ...
    'path_rejections',rejections,'source_and_rng_provenance',true,'manifest_guard',true);
f=fopen(fullfile(out,'checks.json'),'w');assert(f>=0);c=onCleanup(@() fclose(f));
fprintf(f,'%s\n',jsonencode(summary,PrettyPrint=true));
fprintf('PASS production interface: %d trajectories (max %.3g), save/load, wrapper, provenance and safe destinations\n',numel(errors),max(errors));
end

function result=run_wrapper(root,short_cfg)
source=fileread(fullfile(root,'main.m'));
line=regexp(source,'(?m)^cfg=mn3sn_default_config[^\r\n]*','match');assert(numel(line)==1);
source=strrep(source,line{1},'cfg=short_cfg;');
assert(numel(strfind(source,'mn3sn_run(cfg)'))==1);
source=strrep(source,'mn3sn_run(cfg)',"mn3sn_run(cfg,'')");eval(source);
end

function cfg=default_from_root(root)
original=pwd;cleanup=onCleanup(@() cd(original));cd(root);cfg=mn3sn_default_config;
end

function assert_same_trajectory(a,b)
assert(isequal(a.mmx,b.mmx) && isequal(a.mmy,b.mmy) && isequal(a.mmz,b.mmz) && isequal(a.t,b.t));
end

function bytes=read_bytes(file)
f=fopen(file,'rb');assert(f>=0);c=onCleanup(@() fclose(f));bytes=fread(f,Inf,'*uint8');
end

function hash=hash_bytes(bytes)
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function restore(folder,oldpath,oldrng)
cd(folder);path(oldpath);rng(oldrng);
end
