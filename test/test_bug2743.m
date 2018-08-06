function test_bug2743

% MEM 300mb
% WALLTIME 00:10:00

% TEST test_bug2743 ft_mvaranalysis

load(dccnpath('/home/common/matlab/fieldtrip/data/test/test_bug2743.mat'));

cfg = [];
cfg.order = 5;
cfg.toolbox = 'bsmart';
mdata = ft_mvaranalysis(cfg, dataBOWact);

% the issue was a ceil/floor problem in rounding the samples when the tfwin
% had an odd number of samples: this has now been resolved. the purpose of
% the test function is that it should not produce an error in the call to
% ft_mvaranalysis
