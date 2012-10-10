function test_suite = test_bug1168

% TEST test_bug1168
% TEST ft_multiplotTFR

initTestSuite;  % for xUnit


function d = setup

d.cfg = [];
d.cfg.layout = 'CTF151.lay';
d.cfg.lay = ft_prepare_layout(d.cfg);
d.cfg.parameter = 'powspctrm';
d.label = ft_channelselection('MEG', d.cfg.lay.label);
d.nchan = length(d.label);


function test_chan_freq_time(d)
nfreq = 50;
ntime = 100;

freq = [];
freq.label     = d.label;
freq.freq      = 1:nfreq;
freq.time      = 1:ntime;
freq.dimord    = 'chan_freq_time';
freq.powspctrm = randn(d.nchan, nfreq, ntime);

cfg = d.cfg;
ft_multiplotTFR(cfg, freq); % sofar it works


function test_chan_time_freq(d)
nfreq = 50;
ntime = 100;

freq = [];
freq.label = d.label;
freq.freq = 1:nfreq;
freq.time = 1:ntime;
freq.dimord = 'chan_time_freq'; % swap time and freq
freq.powspctrm = randn(d.nchan, ntime, nfreq);

cfg = d.cfg;
ft_multiplotTFR(cfg, freq); % this still works


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following one is more specifically for Linsey
function test_linsey(d)

nchan = length(d.label);
nfreq = 1;
ntime = 100;

freq = [];
freq.label = d.label;
freq.freq = 1:nfreq;
freq.time = 1:ntime;
freq.dimord = 'chan_freq_time';
freq.powspctrm = randn(nchan, nfreq, ntime);

ft_multiplotTFR(d.cfg, freq); % this fails, nothing is plotted, visual
                              % inspection needed


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_unrecognized(d)
% the following rightfully fails because the data structure is not recognized
% as freq structure but if the data structure would have had a freq
% dimension, it should have worked
nax1 = 50;
nax2 = 100;

freq = [];
freq.label = d.label;
freq.ax1 = 1:nax1;
freq.ax2 = 1:nax2;
freq.dimord = 'chan_ax1_ax2';
freq.powspctrm = randn(d.nchan, nax1, nax2);

cfg = d.cfg;
% Test for specific exception using xUnit.
% TODO: find id's for errors (gerp warning -> sort.)
f = @() ft_multiplotTFR(cfg, freq)
assertExceptionThrown(f, 'fieldtrip:dimord')


function test_unknown_dim(d)
% the following rightfully fails because ax2 is not known as dimension
nfreq = 50;
nchan = length(d.label);
nax2 = 100;

freq = [];
freq.label = d.label;
freq.freq = 1:nfreq;
freq.ax2 = 1:nax2;
freq.dimord = 'chan_freq_ax2';
freq.powspctrm = randn(nchan, nfreq, nax2);

% Test for specific exception using xUnit.
% TODO: find id's for errors (gerp warning -> sort.)
f = @() ft_multiplotTFR(cfg, freq)
assertExceptionThrown(f, 'fieldtrip:dimord')
