function test_bug3420

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_statistics_stats ft_timelockstatistics

close all

%%

nchan = 32;
ntrial = 10;
% effect = 0; % in this case the false alarm rate should be close to cfg.alpha, i.e. 5%
effect = 1/15;

timelock = [];
timelock.time = (1:1000)/1000;
timelock.dimord = 'rpt_chan_time';
for i=1:nchan
  timelock.label{i} = num2str(i);
  for j=1:ntrial
    % there is a linear increase over trials
    % the noise varies with the channels
    timelock.trial(j,i,:) = (i/nchan) * randn(1,1000) + j*effect;
  end
end

%%

cfg = [];
cfg.design = 1:ntrial;
cfg.ivar = 1;
cfg.method = 'stats';
cfg.statistic = 'pearson';
pearson = ft_timelockstatistics(cfg, timelock);

figure; imagesc(pearson.stat); colorbar; title('pearson stat');
figure; imagesc(pearson.prob); colorbar; title('pearson prob');

%%

cfg = [];
cfg.design = 1:ntrial;
cfg.ivar = 1;
cfg.method = 'stats';
cfg.statistic = 'kendall';
kendall = ft_timelockstatistics(cfg, timelock);

figure; imagesc(kendall.stat); colorbar; title('kendall stat');
figure; imagesc(kendall.prob); colorbar; title('kendall prob');


%%

cfg = [];
cfg.design = 1:ntrial;
cfg.ivar = 1;
cfg.method = 'stats';
cfg.statistic = 'spearman';
spearman = ft_timelockstatistics(cfg, timelock);

figure; imagesc(spearman.stat); colorbar; title('spearman stat');
figure; imagesc(spearman.prob); colorbar; title('spearman prob');

