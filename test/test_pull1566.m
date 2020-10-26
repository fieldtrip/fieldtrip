function test_pull1566

% MEM 8gb
% WALLTIME 00:20:00
% DEPENDENCY ft_appendfreq append_common ft_selectdata ft_appendtimelock

freq1=[];
freq1.label={'chan1'};
freq1.dimord='chan_freq';
freq1.freq=[1:20];
freq1.powspctrm = rand(1,size(freq1.freq,2));

freq2=freq1;
freq2.label={'chan2'};

cfg=[];
cfg.avgoverfreq='yes';
cfg.keepfreqdim='no';
freq1_s = ft_selectdata(cfg,freq1);
freq2_s = ft_selectdata(cfg,freq2);

cfg=[];
cfg.appenddim = 'chan';
cfg.parameter = 'powspctrm';
freqA = ft_appendfreq(cfg,freq1_s,freq2_s);

% The below assertions are 'incorrect', i.e. they should be throwing an
% error if this issue is solved
%assert(isequal(freqA.freq, [0;0]));
%assert(~isfield(freqA, 'powspctrm'));

% another situation
freq1=[];
freq1.label={'chan1';'chan2'};
freq1.dimord='chan_freq';
freq1.freq=[1:20];
freq1.powspctrm = rand(2,size(freq1.freq,2));

freq2=freq1;
freq2.label={'chan3';'chan4'};

cfg=[];
cfg.avgoverfreq='yes';
cfg.keepfreqdim='no';
freq1_s = ft_selectdata(cfg,freq1);
freq2_s = ft_selectdata(cfg,freq2);

cfg=[];
cfg.appenddim = 'chan';
cfg.parameter = 'powspctrm';
freqB = ft_appendfreq(cfg,freq1_s,freq2_s);

% This scenario yields an OK structure
assert(isfield(freqB, 'powspctrm'));
assert(isequal(freqB.freq, 0));

% another situation
freq1=[];
freq1.label={'chan1'};
freq1.dimord='chan_freq';
freq1.freq=[1:2];
freq1.powspctrm = rand(1,size(freq1.freq,2));

freq2=freq1;
freq2.label={'chan2'};

cfg=[];
cfg.appenddim = 'chan';
cfg.parameter = 'powspctrm';
freqC = ft_appendfreq(cfg,freq1,freq2);

% This scenario seems to correctly append the powerspectra, but the
% freq-field is incorrect, the second assertion should throw an error once
% this has been solved.
assert(isequal(freqC.powspctrm, [freq1.powspctrm;freq2.powspctrm]));
%assert(isequal(freqC.freq, [1 2;1 2]));

% the thing that seems to go wrong is the following:
% line 94 in append_common calls ft_selectdata with multiple inputs (and
% outputs). Once there is just a single channel in the input data, it seems
% that the numeric field freq (and possibly also time, as of yet untested)
% is treated as a 'data field', and appended in the 'unionized' output.

tlck1 = [];
tlck1.label = {'chan1'};
tlck1.time  = [1 2];
tlck1.trial = [0 0];
tlck1.dimord = 'chan_time';

tlck2 = tlck1;
tlck2.label = {'chan2'};
tlck2.trial = [1 1];

cfg = [];
cfg.appenddim = 'chan';
cfg.parameter = 'trial';
tlckA = ft_appendtimelock(cfg, tlck1, tlck2);

% this is incorrect, too:
%assert(isequal(tlckA.time, [1 2;1 2]));

% snippet of code that directly goes into the function in which it goes
% wrong
cfg = [];
cfg.select = 'union';
dataout = cell(1,2);
[dataout{:}] = ft_selectdata(cfg, tlck1, tlck2);
assert(isequal(size(dataout{1}.time), [1 2]));
