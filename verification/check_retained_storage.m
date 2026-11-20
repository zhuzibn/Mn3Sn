function errors=check_retained_storage(before_file)
% Full old captures and stride-one oracles across periodic/open chunk boundaries.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
old=load(before_file,'cases');errors=zeros(1,numel(old.cases));
for k=1:numel(old.cases)
    cfg=old.cases{k}.cfg;reference=old.cases{k}.result;
    actual=mn3sn_run(cfg,'');dense_cfg=cfg;dense_cfg.protocol.savetstep=1;
    dense=mn3sn_run(dense_cfg,'');indices=actual.saved_steps+1;
    for name={'mmx','mmy','mmz'}
        errors(k)=max(errors(k),max(abs(actual.(name{1})-reference.(name{1})),[],'all'));
        assert(isequal(actual.(name{1}),dense.(name{1})(:,:,indices)),...
            'Mn3Sn:Retention','Saving cadence changed the trajectory.');
    end
    assert(errors(k)<=1e-10,'Mn3Sn:Retention','Pre-change trajectory mismatch.');
    for name={'t','saved_steps','actual_updates','pulse_active_updates','config',...
            'requested_config','atomtype','initial','axes','classes'}
        assert(isequal(actual.(name{1}),reference.(name{1})),...
            'Mn3Sn:Retention','Pre-change metadata mismatch: %s',name{1});
    end
end
fprintf('PASS four retention controls: old/new max %.3g, dense subsets exact\n',max(errors));
end
