function check_benchmark(baseline_root,out)
%CHECK_BENCHMARK Comparator mutation tests; no GPU or production rerun.
out=mn3sn_prepare_output(out,baseline_root);
assert(compare_to_baseline(baseline_root,baseline_root),'Baseline self comparison failed.');
cases={'relaxation','driven_nonuniform'};
base_meta=jsondecode(fileread(fullfile(baseline_root,'metadata.json')));
base=cell(1,2);
for k=1:2,s=load(fullfile(baseline_root,[cases{k} '.mat']));base{k}=s.result;end
mutations={'intermediate_rotation','nan','inf','vacancy','shape','mask','axes','classes', ...
    'time','steps','pulse_count','config','initial','schema','summary','missing_field', ...
    'missing_case','missing_inventory','incomplete','mode','metadata_schema'};
for k=1:numel(mutations)
    name=mutations{k};target=fullfile(out,name);mkdir(target);
    r=base{2};m=base_meta;occ=find(r.atomtype==1,1);empty=find(r.atomtype==3,1);
    switch name
        case 'intermediate_rotation'
            [i,j]=ind2sub(size(r.atomtype),occ);angle=1e-4;
            x=r.mmx(i,j,2);y=r.mmy(i,j,2);
            r.mmx(i,j,2)=cos(angle)*x-sin(angle)*y;r.mmy(i,j,2)=sin(angle)*x+cos(angle)*y;
            r.observables=benchmark_validate(r); % Valid norms, final state unchanged.
            assert(isequal(r.mmx(:,:,end),base{2}.mmx(:,:,end)));
        case 'nan',r.mmx(occ)=NaN;
        case 'inf',r.mmy(occ)=Inf;
        case 'vacancy',r.mmz(empty)=1e-6;
        case 'shape',r.mmx=r.mmx(:,:,1:end-1);
        case 'mask',r.atomtype(occ)=3;
        case 'axes',r.axes(occ)=r.axes(occ)+.1;
        case 'classes',r.classes(occ)=mod(r.classes(occ),3)+1;
        case 'time',r.t(2)=r.t(2)+1e-15;
        case 'steps',r.saved_steps(2)=r.saved_steps(2)+1;
        case 'pulse_count',r.pulse_active_updates=r.pulse_active_updates+1;
        case 'config',r.config.alp=r.config.alp*2;
        case 'initial',r.initial(occ)=r.initial(occ)+.1;
        case 'schema',r.schema_version=999;
        case 'summary',r.observables.mean(2,1)=r.observables.mean(2,1)+1e-4;
        case 'missing_field',r=rmfield(r,'mmx');
        case 'missing_inventory',m.cases=rmfield(m.cases,'driven_nonuniform');
        case 'incomplete',m.status='running';
        case 'mode',m.mode='wrong';
        case 'metadata_schema',m.schema_version=999;
    end
    result=base{1};save(fullfile(target,'relaxation.mat'),'result');
    if ~strcmp(name,'missing_case'),result=r;save(fullfile(target,'driven_nonuniform.mat'),'result');end
    write_json(fullfile(target,'metadata.json'),m);
    [passed,report]=compare_to_baseline(target,baseline_root);
    assert(~passed,'Mutation was accepted: %s',name);
    if strcmp(name,'intermediate_rotation')
        assert(isfield(report,'cases') && isfield(report.cases,'driven_nonuniform') && ...
            report.cases.driven_nonuniform.trajectory_max_abs>1e-10, ...
            'Rotation did not exercise trajectory comparison.');
    end
end
assert(~compare_to_baseline(baseline_root,fullfile(out,'nan')),'Invalid baseline accepted.');
% A changed source revision is allowed; bad timing suppresses speedup only.
target=fullfile(out,'provenance_and_timing');mkdir(target);m=base_meta;
m.source.test_revision='future-change';
for k=1:2
    result=base{k};
    if k==1,result.timing.warm_seconds=[NaN -1 0];
    else,result.timing.median_seconds=result.timing.median_seconds+1i;end
    save(fullfile(target,[cases{k} '.mat']),'result');
end
write_json(fullfile(target,'metadata.json'),m);
[passed,report]=compare_to_baseline(target,baseline_root);
assert(passed && report.source_changed && ~report.cases.relaxation.performance_comparable && ...
    ~report.cases.driven_nonuniform.performance_comparable);
% Existing output must be rejected without touching any files.
old=fileread(fullfile(baseline_root,'metadata.json'));rejected=false;
try,run_baseline_benchmark(base_meta.mode,baseline_root);catch,rejected=true;end
assert(rejected && strcmp(old,fileread(fullfile(baseline_root,'metadata.json'))));
root=fileparts(fileparts(mfilename('fullpath')));rejected=false;
try,run_baseline_benchmark(base_meta.mode,fullfile(root,'benchmark','forbidden-output'));catch,rejected=true;end
assert(rejected && ~isfolder(fullfile(root,'benchmark','forbidden-output')));
% New destinations within either registered or explicit baseline must fail.
for target={fullfile(baseline_root,'new-candidate'), ...
        fullfile(baseline_root,'unused','..','alias-candidate')}
    rejected=false;
    try,run_baseline_benchmark(base_meta.mode,target{1},baseline_root);
    catch err,rejected=strcmp(err.identifier,'Mn3Sn:ProtectedBaseline');end
    assert(rejected && ~isfolder(target{1}));
end
missing=fullfile(out,'missing-input');target=fullfile(out,'missing-input-output');rejected=false;
try,run_baseline_benchmark(base_meta.mode,target,missing);catch,rejected=true;end
assert(rejected && ~isfolder(target));
% Missing-manifest failure is exercised in a copy, never by moving the reference.
copy=fullfile(out,'missing-manifest-source');mkdir(copy);mkdir(fullfile(copy,'benchmark'));
copyfile(fullfile(root,'benchmark','run_baseline_benchmark.m'),fullfile(copy,'benchmark'));
oldpath=path;oldpwd=pwd;context=onCleanup(@() restore_context(oldpwd,oldpath));
cd(copy);addpath(fullfile(copy,'benchmark'),'-begin');
target=fullfile(out,'missing-manifest-output');rejected=false;
try,run_baseline_benchmark(base_meta.mode,target);
catch err,rejected=strcmp(err.identifier,'Mn3Sn:MissingBaseline');end
assert(rejected && ~isfolder(target));
cd(oldpwd);path(oldpath);clear context;
files=dir(fullfile(root,'benchmark','*.m'));issues=struct;
for k=1:numel(files)
    messages=checkcode(fullfile(files(k).folder,files(k).name),'-id');
    issues.(matlab.lang.makeValidName(files(k).name))=messages;
    assert(~any(strcmp({messages.id},'SYNER')),'MATLAB syntax error: %s',files(k).name);
end
summary=struct('passed',true,'mutation_rejections',numel(mutations), ...
    'invalid_baseline_rejected',true,'source_change_allowed',true, ...
    'invalid_timing_excluded',true,'path_guards_passed',true, ...
    'missing_references_rejected',true,'baseline_subtrees_protected',true,'checkcode',issues);
write_json(fullfile(out,'checks.json'),summary);
fprintf('PASS benchmark checks: %d mutations, invalid baseline, timing/source policy, safe paths.\n',numel(mutations));
end

function restore_context(folder,oldpath)
cd(folder);path(oldpath);
end

function write_json(file,value)
f=fopen(file,'w');assert(f>=0);c=onCleanup(@() fclose(f));fprintf(f,'%s\n',jsonencode(value,PrettyPrint=true));
end
