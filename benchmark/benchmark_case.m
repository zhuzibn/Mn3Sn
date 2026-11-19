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
assert(strcmpi(char(java.io.File(fileparts(which('mn3sn_initialize'))).getCanonicalPath()), ...
    char(java.io.File(root).getCanonicalPath())),'Initialization resolved outside the selected source tree.');
cfg=mn3sn_default_config;cfg.protocol=selected_protocol;
state=mn3sn_initialize(cfg);c=state.config;
assert(c.rk4==1 && c.bc==0 && ~c.thermalenable && ~c.dipolemode && ~c.enablefixedge);
assert(~c.loadstartm && ~c.load_fixed_atom_distrib && ~c.save_fixed_atom_distrib && ~c.dwcalc);
assert(c.SOT_DLT==1 && c.SOT_FLT==0 && c.STT_DLT==0 && c.STT_FLT==0);
result=struct('schema_version',1,'mode',mode,'case_name',name);
result.atomtype=state.atomtype;result.axes=state.axes;
result.generated_initial=state.generated_initial;result.config=c;
result.classes=state.classes;
occ=find(state.atomtype==1);N=numel(occ);
assert(N==c.natomW*c.natomL*3/4 && all(ismember(state.atomtype(:),[1 3])));
assert(all(arrayfun(@(q)nnz(state.classes==q)==N/3,1:3)));
angles=(c.initial_theta_deg+(0:2)*120)*pi/180;
my_init=state.generated_initial(:,:,2);mz_init=state.generated_initial(:,:,3);
for q=1:3
    mask=state.classes==q;
    assert(max(abs(my_init(mask)-sin(angles(q))))<1e-14 && ...
        max(abs(mz_init(mask)-cos(angles(q))))<1e-14);
end
if ~isempty(frozen)
    for key={'schema_version','mode','case_name','atomtype','axes','classes','generated_initial','config'}
        assert(isequaln(result.(key{1}),frozen.(key{1})), ...
            'Frozen input configuration mismatch: %s',key{1});
    end
end
mx_init=state.initial(:,:,1);my_init=state.initial(:,:,2);mz_init=state.initial(:,:,3);
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
cfg.initial_state=result.initial;
production=mn3sn_run(cfg,'');
assert(isequaln(production.config,result.config) && isequal(production.initial,result.initial));
result.setup_seconds=result.setup_seconds+production.setup_seconds;
for key={'integration_seconds','mmx','mmy','mmz','saved_steps','t','actual_updates','pulse_active_updates'}
    result.(key{1})=production.(key{1});
end
result.observables=benchmark_validate(result);
end
