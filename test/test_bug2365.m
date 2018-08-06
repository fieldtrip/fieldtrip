function test_bug2365

% MEM 500mb
% WALLTIME 00:10:00

% TEST ft_freqanalysis

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2365.mat');
load(filename); % loads variable 'data'

cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 1:40;
cfg.t_ftimwin = 4 ./ cfg.foi;
cfg.toi = 'all';
cfg.keeptrials = 'yes';
cfg.output = 'pow';

convol = ft_freqanalysis(cfg, data);

cfg.output = 'fourier'; % does not work
convol2 = ft_freqanalysis(cfg, data);
