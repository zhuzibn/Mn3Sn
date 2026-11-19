%% Mn3Sn deterministic GPU simulation; see README for inputs and result fields.
cfg=mn3sn_default_config('relaxation'); % e.g. ('fig4a'), ('fig5a',300)
result=mn3sn_run(cfg); % fresh external run directory; use mn3sn_run(cfg,path) to choose
fprintf('Saved Mn3Sn result: %s\n',fullfile(result.output_root,'result.mat'));
