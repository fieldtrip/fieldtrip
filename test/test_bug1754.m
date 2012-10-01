load test_bug1754

cfg = [];
cfg.baseline = [-0.3 -0.2];
cfg.baselinetype = 'relchange';
freqb = ft_freqbaseline(cfg, freq);

figure;imagesc(squeeze(freqb.powspctrm(1,:,:)));
