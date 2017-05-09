function test_ft_appendfreq

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_appendfreq

% make some dummy frequency structures
freq1.label = {'1';'2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(2,10,5);

cfg = [];
cfg.parameter = 'powspctrm';

freq2         = freq1;
cfg.appenddim = 'rpt';
freqrpt       = ft_appendfreq(cfg, freq1, freq2);
cfg.appenddim = 'auto';
freqrptauto   = ft_appendfreq(cfg, freq1, freq2);

freq2         = freq1;
freq2.label   = {'3';'4'};
cfg.appenddim = 'chan';
freqchan      = ft_appendfreq(cfg, freq1, freq2);
cfg.appenddim = 'auto';
freqchanauto  = ft_appendfreq(cfg, freq1, freq2);

freq2         = freq1;
freq2.freq    = 11:20;
cfg.appenddim = 'freq';
freqfreq      = ft_appendfreq(cfg, freq1, freq2);
if ~isfield(freqfreq,'freq')
  error('freq field is not appeneded: see bugs 1984 and 2187');
end

cfg.appenddim = 'auto';
freqfreqauto  = ft_appendfreq(cfg, freq1, freq2);
if ~isfield(freqfreqauto,'freq')
  error('freq field is not appeneded: see bugs 1984 and 2187');
end

freq2         = freq1;
freq2.time    = 6:10;
cfg.appenddim = 'time';
freqtime      = ft_appendfreq(cfg, freq1, freq2);
if ~isfield(freqfreqauto,'freq')
  error('freq field is not appeneded: see bugs 1984 and 2187');
end

cfg.appenddim = 'auto';
freqtimeauto  = ft_appendfreq(cfg, freq1, freq2);
if ~isfield(freqtimeauto,'freq')
  error('freq field is not appeneded: see bugs 1984 and 2187');
end

% now test for numerical inaccurracies, should concatenate across 'rpt'
freq2          = freq1;
freq2.time     = freq1.time+0.0000001;
cfg.appenddim  = 'auto';
freqrpt2       = ft_appendfreq(cfg, freq1, freq2);
if ~isfield(freqrpt2,'freq')
  error('freq field is not appeneded: see bugs 1984 and 2187');
end

%% test for data with labels shuffled around
freq1 = [];
freq1.label = {'1';'2'};
freq1.freq  = 1:10;
freq1.time  = 1:5;
freq1.dimord = 'chan_freq_time';
freq1.powspctrm = randn(2,10,5);

freq2 = freq1;
freq2.label = {'2';'1'};
freq2.powspctrm = randn(2,10,5);

freqshuffled = ft_appendfreq(cfg, freq1, freq2);

if ~strcmp(freqshuffled.dimord, 'rpt_chan_freq_time')
  error('unexpected dimord when appending freqs with differently permuted labels');
end
if ~isfield(freqshuffled,'freq')
  error('freq field is not appeneded: see bugs 1984 and 2187');
end

% now check whether channels were correctly appended
[a,b] = match_str(freq1.label, freqshuffled.label);
x1 = freq1.powspctrm(a(1),:,:);
y1 = freqshuffled.powspctrm(1,b(1),:,:);

[a,b] = match_str(freq2.label, freqshuffled.label);
x2 = freq2.powspctrm(a(1),:,:);
y2 = freqshuffled.powspctrm(2,b(1),:,:);

if ~all(x1(:) == y1(:)) || ~all(x2(:) == y2(:))
  error('data was wrongly appended when channel labels are differently ordered in input arguments');
end

