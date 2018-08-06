function test_bug3182

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqanalysis ft_freqdescriptives

%% generate some data

data          = [];
for i=1:10
  data.trial{i} = randn(2,1000);
  data.time{i}  = (1:1000)/1000;
end
data.label    = {'chan1';'chan2'};
data.fsample  = 1000;
data.cfg.trl  = zeros(2,3);

%% compute

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 4;
cfg.foilim    = [8 40];

cfg.keeptrials = 'no';
freq1 = ft_freqanalysis(cfg, data);

cfg.keeptrials = 'yes';
freq2 = ft_freqanalysis(cfg, data);

cfg = [];
freq2fd = ft_freqdescriptives(cfg, freq2);

%% compare

if max(abs(freq1.powspctrm(:)-freq2fd.powspctrm(:))./norm(freq1.powspctrm(:)+freq2fd.powspctrm(:)))>1e-10
  % small errors are expected due to accumulation of numerical errors
  error('they are different');
end


