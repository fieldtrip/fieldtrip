function test_bug1508

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqanalysis

% Stan reported a strange error caused by the following:
% if the cfg.channel in a call to ft_freqanalysis contains channels that
% are not in the data, the function crashes on the second trial

%% Try to reproduce
data = [];
data.trial = {randn(3,100) randn(3,100)};
data.time  = {1:100 1:100};
data.label = {'chan1';'chan2';'chan3'};

cfg = [];
cfg.method  = 'mtmfft';
cfg.channel = {'chan4'};
cfg.output  = 'pow';
cfg.taper   = 'hanning';

try
  freq = ft_freqanalysis(cfg, data);
catch err
  if strcmp(err.message,'no channels were selected')
    % error is expected handling 
  else
    error('test doesn''t result in expected error message')
  end
end
