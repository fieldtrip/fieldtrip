function test_bug2086

% MEM 1500mb
% WALLTIME 00:03:01

% TEST test_bug2086
% TEST ft_databrowser

warning('this test should not run automatically');
return

load(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2086.mat'));

cfg = [];
cfg.continuous = 'no';
cfg.preproc.channel = {'all'  };
cfg.preproc.demean = 'yes';
cfg.viewmode = 'butterfly';
cfg.ylim = [-200 200];
cfg.channel = 'all';
cfg = ft_databrowser(cfg,data);
