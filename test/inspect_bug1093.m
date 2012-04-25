function cfgnew = test_bug1093

% TEST test_bug1093 ft_artifact_zvalue

% the original bug was that ft_artifact_zvalue does not adjust the detected
% artifacts when the threshold is adjusted

load test_bug1027; % use the data for bug1027

cfg = [];
cfg.artfctdef.zvalue.channel = {'EOG061'};
cfg.artfctdef.zvalue.cutoff  = 4;
cfg.artfctdef.zvalue.interactive = 'yes';
cfgnew = ft_artifact_zvalue(cfg, data);