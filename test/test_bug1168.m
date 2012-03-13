function test_bug1168

% TEST test_bug1168
% TEST ft_multiplotTFR

cfg = [];
cfg.layout = 'CTF151.lay';
lay = ft_prepare_layout(cfg);

label = ft_channelselection('MEG', lay.label);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nchan = length(label);
nfreq = 50;
ntime = 100;

freq = [];
freq.label = label;
freq.freq = 1:nfreq;
freq.time = 1:ntime;
freq.dimord = 'chan_freq_time';
freq.powspctrm = randn(nchan, nfreq, ntime);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = lay;
ft_multiplotTFR(cfg, freq) % sofar it works

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = [];
freq.label = label;
freq.freq = 1:nfreq;
freq.time = 1:ntime;
freq.dimord = 'chan_time_freq'; % swap time and freq
freq.powspctrm = randn(nchan, ntime, nfreq);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.layout = lay;
ft_multiplotTFR(cfg, freq) % this still works


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following one is more specifically for Linsey
failed = true;

nchan = length(label);
nfreq = 1;
ntime = 100;

freq = [];
freq.label = label;
freq.freq = 1:nfreq;
freq.time = 1:ntime;
freq.dimord = 'chan_freq_time';
freq.powspctrm = randn(nchan, nfreq, ntime);
try
  ft_multiplotTFR(cfg, freq) % this fails, nothing is plotted, visual inspection needed
  failed = false;
catch
  failed = true;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


try
  % the following rightfully fails because the data structure is not recognized as freq structure
  % but if the data structure would have had a freq dimension, it should have worked
  nax1 = 50;
  nax2 = 100;
  
  freq = [];
  freq.label = label;
  freq.ax1 = 1:nax1;
  freq.ax2 = 1:nax2;
  freq.dimord = 'chan_ax1_ax2';
  freq.powspctrm = randn(nchan, nax1, nax2);
  
  cfg = [];
  cfg.parameter = 'powspctrm';
  cfg.layout = lay;
  ft_multiplotTFR(cfg, freq)
  failed = false;
catch
  failed = failed & true;
end

try
  % the following rightfully fails because ax2 is not known as dimension
  nfreq = 50;
  nax2 = 100;
  
  freq = [];
  freq.label = label;
  freq.freq = 1:nfreq;
  freq.ax2 = 1:nax2;
  freq.dimord = 'chan_freq_ax2';
  freq.powspctrm = randn(nchan, nfreq, nax2);
  
  cfg = [];
  cfg.parameter = 'powspctrm';
  cfg.layout = lay;
  ft_multiplotTFR(cfg, freq)
  failed = false;
catch
  failed = failed & true;
end

if ~failed
  error('either of the two cases which should fail works')
end
  
