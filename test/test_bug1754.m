function test_bug1754
% Note: new tests were added when bug #1754 was fixed. This script is now
% reduntant, and can be removed.
load test_bug1754

cfg = [];
cfg.baseline = [-0.3 -0.2];
cfg.baselinetype = 'relchange';
freqb = ft_freqbaseline(cfg, freq);

figure;imagesc(squeeze(freqb.powspctrm(1,:,:)));
