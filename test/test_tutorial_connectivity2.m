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
  dat = (dat-repmat(mean(dat,2),[1 nSamples]))./repmat(std(dat,[],2),[1 nSamples]);
  data.trial{k} = mixing * dat;
  data.time{k}  = (0:nSamples-1)./fsample;
end
data.label = {'chan1' 'chan2'}';

figure;plot(dat'+repmat([0 1 2],[nSamples 1]));
title('original ''sources''');

figure;plot((mixing*dat)'+repmat([0 1],[nSamples 1])); 
axis([0 1000 -1 2]);
set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
title('mixed ''sources''');

%% do spectral analysis
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [0 200];
cfg.tapsmofrq = 5;
freq = ft_freqanalysis(cfg, data);
fd   = ft_freqdescriptives([], freq);

figure;plot(fd.freq, fd.powspctrm);
set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
title('powerpectrum');

%% compute connectivity
cfg = [];
cfg.method = 'granger';
g = ft_connectivityanalysis(cfg, freq);
cfg.method = 'coh';
c = ft_connectivityanalysis(cfg, freq);

%% visualize
cfg = [];
cfg.parameter = 'grangerspctrm';
ft_connectivityplot(cfg, g);
title('granger causality');
cfg.parameter = 'cohspctrm';
ft_connectivityplot(cfg, c);
title('coherence');

