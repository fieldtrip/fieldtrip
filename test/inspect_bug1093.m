function inspect_bug1093

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_artifact_zvalue
% DATA private

% the original bug was that ft_artifact_zvalue does not adjust the detected
% artifacts when the threshold is adjusted

% use the data for bug1027
cd(dccnpath('/project/3031000.02/test'))
load bug1027.mat

cfg = [];
cfg.artfctdef.zvalue.channel = {'EOG061'};
cfg.artfctdef.zvalue.cutoff  = 4;
cfg.artfctdef.zvalue.interactive = 'yes';
cfgnew = ft_artifact_zvalue(cfg, data);
