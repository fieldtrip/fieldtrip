function cfgnew = test_bug1093

% TEST test_bug1093 ft_artifact_zvalue

% the original bug was that ft_artifact_zvalue does not adjust the detected
% artifacts when the threshold is adjusted

% use the data for bug1027
cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1027.mat

cfg = [];
cfg.artfctdef.zvalue.channel = {'EOG061'};
cfg.artfctdef.zvalue.cutoff  = 4;
cfg.artfctdef.zvalue.interactive = 'yes';
cfgnew = ft_artifact_zvalue(cfg, data);
