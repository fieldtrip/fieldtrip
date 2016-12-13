function test_ft_freqbaseline

% MEM 200mb
% WALLTIME 00:10:00

% TEST test_ft_freqbaseline 
% TEST ft_freqbaseline

% generate some data
freq = [];
freq.time = (-460:10:1000)./1000;
freq.freq = 10:2:100;
freq.label = {'chan01';'chan02'};
freq.dimord = 'chan_freq_time';
freq.powspctrm = rand(2, numel(freq.freq), numel(freq.time))./5;
nfreq = numel(freq.freq);
for k = 1:nfreq
  indx = 47+((-nfreq+k):(nfreq-k));
  freq.powspctrm(:,k,indx) = freq.powspctrm(:,k,indx)+1;
end

cfg = [];
cfg.baseline = [-inf -0.1];
cfg.baselinetype = 'absolute';
fb1 = ft_freqbaseline(cfg, freq);

dtime = 0.01;
cfg.baseline = [ones(numel(freq.freq),1).*-inf -((numel(freq.freq)):-1:1)'.*dtime];
fb2 = ft_freqbaseline(cfg, freq);

