function test_bug2473

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_databrowser

load(dccnpath('/home/common/matlab/fieldtrip/data/test/avgFIC.mat'));

cfg = [];
cfg.preproc.bpfilter = 'yes';
cfg.preproc.bpfreq = [3 45];
ft_databrowser(cfg, avgFIC);

% Opening the preproc window in the databrowser and clicking on save and
% close produces an error because the cfg.preproc.bpfreq value is not on
% the same line anymore.

end
