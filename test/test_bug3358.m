function test_bug3358

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY

d.cfg = [];
d.cfg.layout = 'CTF151.lay';
d.cfg.lay = ft_prepare_layout(d.cfg);
d.cfg.parameter = 'powspctrm';
d.label = ft_channelselection('MEG', d.cfg.lay.label);
d.nchan = length(d.label);

nchan = d.nchan;
nfreq = 50;
ntime = 100;

%%

freq = [];
freq.label     = d.label;
freq.freq      = 1:nfreq;
freq.time      = 1:ntime;
freq.dimord    = 'chan_freq_time';
freq.powspctrm = randn(nchan, nfreq, ntime);

cfg = [];
cfg.layout = 'CTF151';
figure
ft_multiplotTFR(cfg, freq);

%%

freq = [];
freq.label     = d.label;
freq.freq      = 1:nfreq;
freq.time      = 1:ntime;
freq.dimord    = 'chan_time_freq';
freq.powspctrm = randn(nchan, ntime, nfreq);

cfg = [];
cfg.layout = 'CTF151';
figure
ft_multiplotTFR(cfg, freq);


%%
% This one is too weird

% freq = [];
% freq.label     = d.label;
% freq.freq      = 1:nfreq;
% freq.time      = 1:ntime;
% freq.dimord    = 'time_freq_chan';
% freq.powspctrm = randn(ntime, nfreq, nchan);
%
% cfg = [];
% cfg.layout = 'CTF151';
% ft_multiplotTFR(cfg, freq);

%%
% This one is too weird

% freq = [];
% freq.label     = d.label;
% freq.aa        = 1:nfreq;
% freq.bb        = 1:ntime;
% freq.dimord    = 'chan_aa_bb';
% freq.powspctrm = randn(nchan, ntime, nfreq);
%
% cfg = [];
% cfg.layout = 'CTF151';
% ft_multiplotTFR(cfg, freq);
