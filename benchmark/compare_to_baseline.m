function [passed,report]=compare_to_baseline(candidate_root,baseline_root)
%COMPARE_TO_BASELINE Read-only full-trajectory regression and timing comparison.
% [passed,report]=compare_to_baseline(candidate_directory,baseline_directory)
% assert(passed) in MATLAB -batch to fail automation on any correctness error.
report=struct('schema_version',1,'passed',false,'candidate',char(candidate_root), ...
    'baseline',char(baseline_root),'tolerance',1e-10,'error','');
passed=false;
try
    cm=jsondecode(fileread(fullfile(candidate_root,'metadata.json')));
    bm=jsondecode(fileread(fullfile(baseline_root,'metadata.json')));
    assert(cm.schema_version==1 && bm.schema_version==1,'Unsupported suite schema.');
    assert(strcmp(cm.status,'complete') && strcmp(bm.status,'complete'),'Suite is incomplete.');
    assert(strcmp(cm.mode,bm.mode) && ismember(cm.mode,{'quick','current'}),'Suite mode mismatch.');
    report.mode=cm.mode;
    report.same_environment=isequaln(cm.gpu,bm.gpu) && strcmp(cm.matlab,bm.matlab) && ...
        strcmp(cm.os,bm.os) && isequaln(cm.toolboxes,bm.toolboxes);
    report.source_changed=~isequaln(cm.source,bm.source);
    cases={'relaxation','driven_nonuniform'};
    assert(isequal(sort(fieldnames(cm.cases)),sort(cases')) && ...
        isequal(sort(fieldnames(bm.cases)),sort(cases')),'Incomplete case inventory.');
    for k=1:2
        name=cases{k};a=load(fullfile(candidate_root,[name '.mat']),'result');
        b=load(fullfile(baseline_root,[name '.mat']),'result');a=a.result;b=b.result;
        assert(strcmp(a.case_name,name) && strcmp(b.case_name,name) && ...
            strcmp(a.mode,cm.mode) && strcmp(b.mode,bm.mode),'Case identity mismatch.');
        ao=benchmark_validate(a);bo=benchmark_validate(b);
        for key={'schema_version','mode','case_name','config','atomtype','axes', ...
                'classes','generated_initial','initial','saved_steps','actual_updates','pulse_active_updates'}
            assert(isequaln(a.(key{1}),b.(key{1})),'Incompatible %s: %s',name,key{1});
        end
        check_observables(a.observables,ao);check_observables(b.observables,bo);
        s=struct;
        s.trajectory_max_abs=max([max(abs(a.mmx-b.mmx),[],'all'), ...
            max(abs(a.mmy-b.mmy),[],'all'),max(abs(a.mmz-b.mmz),[],'all')]);
        s.mean_max_abs=max(abs(ao.mean-bo.mean),[],'all');
        s.class_mean_max_abs=max(abs(ao.class_mean-bo.class_mean),[],'all');
        s.time_max_abs=max(abs(a.t-b.t));
        s.norm_error_candidate=ao.norm_error;s.norm_error_baseline=bo.norm_error;
        s.passed=s.trajectory_max_abs<=report.tolerance && s.mean_max_abs<=report.tolerance && ...
            s.class_mean_max_abs<=report.tolerance && s.time_max_abs<=8*eps(a.config.runtime);
        [at,av]=timing(a.timing);[bt,bv]=timing(b.timing);
        s.candidate_median_seconds=at;s.baseline_median_seconds=bt;
        s.performance_comparable=report.same_environment && av && bv;
        s.speedup=NaN;if s.performance_comparable,s.speedup=bt/at;end
        report.cases.(name)=s;
        fprintf('%s: trajectory %.3g, mean %.3g, classes %.3g, PASS %d\n', ...
            name,s.trajectory_max_abs,s.mean_max_abs,s.class_mean_max_abs,s.passed);
        if s.performance_comparable
            fprintf('  warm medians baseline %.6f s / candidate %.6f s; speedup %.3fx (informational)\n',bt,at,s.speedup);
        else
            fprintf('  timing not comparable: different environment or invalid timing metadata\n');
        end
    end
    passed=report.cases.relaxation.passed && report.cases.driven_nonuniform.passed;
    if ~passed,report.error='Trajectory or observable tolerance exceeded.';end
catch err
    report.error=err.message;fprintf('FAIL: %s\n',err.message);
end
report.passed=passed;fprintf('Benchmark PASS: %d\n',passed);
end

function check_observables(stored,recomputed)
for key={'mean','class_mean','norm_error','vacancy_max'}
    x=stored.(key{1});y=recomputed.(key{1});
    assert(isnumeric(x) && isreal(x) && isequal(size(x),size(y)) && ...
        all(isfinite(x),'all') && max(abs(x-y),[],'all')<=1e-14, ...
        'Stored observable inconsistent with arrays: %s',key{1});
end
end

function [seconds,valid]=timing(t)
seconds=NaN;valid=false;
if ~isstruct(t) || ~isfield(t,'warm_seconds') || ~isfield(t,'median_seconds'),return;end
v=t.warm_seconds;
if ~isnumeric(v) || ~isreal(v) || numel(v)~=3 || any(~isfinite(v)) || any(v<=0),return;end
seconds=median(v);
valid=isnumeric(t.median_seconds) && isreal(t.median_seconds) && ...
    isscalar(t.median_seconds) && isfinite(t.median_seconds) && ...
    abs(seconds-t.median_seconds)<=8*eps(seconds);
if ~valid,seconds=NaN;end
end
