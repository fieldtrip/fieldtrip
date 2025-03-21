function test_bug1754

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqbaseline
% DATA private

% Note: new tests were added when bug #1754 was fixed. This script is now
% reduntant, and can be removed.

cd(dccnpath('/project/3031000.02/test'))
load bug1754.mat

cfg = [];
cfg.baseline = [-0.3 -0.2];
cfg.baselinetype = 'relchange';
freqb = ft_freqbaseline(cfg, freq);

figure
imagesc(squeeze(freqb.powspctrm(1,:,:)));

