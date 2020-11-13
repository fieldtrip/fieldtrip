function test_bug103

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_singleplotER

freq.freq       = 1:1:100;
freq.powspctrm  = randn(size(freq.freq)).^2;
freq.label      = {'chan1'};
freq.dimord     = 'chan_freq';

cfg = [];
figure; ft_singleplotER(cfg, freq);

% save to a temporary file
filename = [tempname,'.mat'];
save(filename, 'freq');

try
  cfg = [];
  cfg.inputfile = filename;
  figure; ft_singleplotER(cfg);
  delete(filename);
catch ME
  delete(filename);
  rethrow(ME);
end
