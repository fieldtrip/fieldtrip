function test_pull1566

% MEM 8gb
% WALLTIME 00:20:00
% DEPENDENCY ft_appendfreq append_common ft_selectdata

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
assert(isequal(freqA.freq, [0;0]));
assert(~isfield(freqA, 'powspctrm'));

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
assert(isequal(freqC.freq, [1 2;1 2]));
