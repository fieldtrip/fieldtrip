function test_bug1027

% TEST test_bug1027 ft_artifact_zvalue

% the original bug was that ft_artifact_zvalue could not deal with variable
% length trials

load test_bug1027

cfg = [];
cfg.artfctdef.zvalue.channel = {'EOG061'};
cfg.artfctdef.zvalue.cutoff  = 4;
cfg.artfctdef.zvalue.interactive = 'yes';
ft_artifact_zvalue(cfg, data);