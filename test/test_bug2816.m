function test_bug2816

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_selectdata ft_timelockanalysis getdimord getdimsiz

data.label = {'TETFP09'  'TETFP10'  'TETFP11'  'TETFP12'};
for i=1:142
  data.time{i} = (1:1000)./1000;
  data.trial{i} = randn(4,1000);
end
data.fsample = 1000;
data.sampleinfo = rand(142,2);
data.trialinfo = rand(142,1);
% data.hdr is not needed
% data.cfg is not needed

cfg = [];
cfg.covariancewindow = [0.1 0.2];
cfg.keeptrials = 'yes';
cfg.removemean = 'yes';
cfg.covariance = 'yes';
cfg.channel = 'TETFP10';
cfg.trials = [1:2:130];

avg = ft_timelockanalysis(cfg, data);
