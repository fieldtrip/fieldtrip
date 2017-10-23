function test_bug1754

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqbaseline

% Note: new tests were added when bug #1754 was fixed. This script is now
% reduntant, and can be removed.

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1754.mat

cfg = [];
cfg.baseline = [-0.3 -0.2];
cfg.baselinetype = 'relchange';
freqb = ft_freqbaseline(cfg, freq);

figure
imagesc(squeeze(freqb.powspctrm(1,:,:)));

