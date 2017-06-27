function test_bug1027

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_artifact_zvalue

% the original bug was that ft_artifact_zvalue could not deal with variable
% length trials

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'));
load bug1027.mat

cfg = [];
cfg.artfctdef.zvalue.channel = {'EOG061'};
cfg.artfctdef.zvalue.cutoff  = 4;
cfg.artfctdef.zvalue.interactive = 'no';
ft_artifact_zvalue(cfg, data);
