function test_bug2137

% MEM 1gb
% WALLTIME 00:30:00
% DEPENDENCY
% DATA public


% the following lines are for interactive/manual testing
if false
  restoredefaultpath
  clear all
  addpath(dccnpath('/home/common/matlab/fieldtrip'));
  ft_defaults
end

cfg = [];
cfg.dataset = dccnpath('/project/3031000.02/external/download/test/ctf/Subject01.ds');
data = ft_preprocessing(cfg);

cfg            = [];
cfg.detrend    = 'no';
cfg.resamplefs = 150;
data           = ft_resampledata(cfg, data);

cfg            = [];
cfg.channel    = {'MEG'};
cfg.numcomponent = 60;
comp           = ft_componentanalysis(cfg, data);

