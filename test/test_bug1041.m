function test_bug1041

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqdescriptives memtic memtoc

% this is a bug that Jorn reported but that I am not able to reproduce
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1041

% generate some data
data          = [];
data.trial{1} = randn(2,1000);
data.time{1}  = [0:0.001:0.999] - 0.5;
data.trial{2} = randn(2,1000);
data.time{2}  = [0.001:0.001:1] - 0.5;
data.label    = {'chan1';'chan2'};
data.fsample  = 1000;
data.cfg.trl  = zeros(2,3);

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 4;
cfg.foilim    = [8 40];
freq = ft_freqanalysis(cfg, data);

cfg = [];
fd = ft_freqdescriptives(cfg, freq);

if isnan(fd.cfg.callinfo.procmem)
  error('memtoc failed for ft_freqdescriptives');
end

