function test_tutorial_connectivity2

% contains the code that is in the second part of the ocnnectivity
% tutorial. Purpose: simulate linearly mixed data, and show that dependent
% on the parameter settings the connectivity can go wild.

%% create some instantaneously mixed data
nTrials  = 100;
nSamples = 1000;
fsample  = 1000;
mixing   = [0.8 0.2 0;0 0.2 0.8]; % the middle signal will be common pick up


data = [];
data.trial = cell(1,nTrials);
data.time  = cell(1,nTrials);
for k = 1:nTrials
  dat = randn(3, nSamples);
  dat(2,:) = ft_preproc_bandpassfilter(dat(2,:), 1000, [15 25]);
  dat([1 3],:) = dat([1 3],:)./5; % do some scaling (does not matter what)
  data.trial{k} = mixing * dat;
  data.time{k}  = (0:nSamples-1)./fsample;
end
data.label = {'chan1' 'chan2'}';

%% do spectral analysis
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [0 200];
cfg.tapsmofrq = 5;
freq = ft_freqanalysis(cfg, data);

%% compute connectivity
cfg = [];
cfg.method = 'granger';
g = ft_connectivityanalysis(cfg, freq);
cfg.method = 'coh';
c = ft_connectivityanalysis(cfg, freq);

%% visualize
cfg = [];
cfg.parameter = 'grangerspctrm';
figure;ft_connectivityplot(cfg, g);
cfg.parameter = 'cohspctrm';
figure;ft_connectivityplot(cfg, c);
