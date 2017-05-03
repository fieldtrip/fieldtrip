function test_bug3267

% WALLTIME 00:10:00
% MEM 3gb

% this script tests that the trialinfo is properly dealt with in the ft_appendxxx functions

dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];                                           % empty configuration
cfg.dataset                 = dataset;              % name of CTF dataset
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 3;                    % trigger value for fully incongruent (FIC)
cfg = ft_definetrial(cfg);
dataFIC = ft_preprocessing(cfg);

cfg = [];                                           % empty configuration
cfg.dataset                 = dataset;              % name of CTF dataset
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 9;                    % trigger value for fully congruent (FC)
cfg = ft_definetrial(cfg);
dataFC = ft_preprocessing(cfg);

cfg = [];                                           % empty configuration
cfg.dataset                 = dataset;              % name of CTF dataset
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = [3 9];
cfg = ft_definetrial(cfg);
dataAll = ft_preprocessing(cfg);

cfg = [];
dataAppend = ft_appenddata(cfg, dataFIC, dataFC);

%% compare the two versions

% the order should be different
assert(~isequal(dataAll.trialinfo,  dataAppend.trialinfo));
assert(~isequal(dataAll.sampleinfo, dataAppend.sampleinfo));

% sort them according to the trial order in the original experiment
[dum, order1] = sort(dataAll.sampleinfo(:,1));
[dum, order2] = sort(dataAppend.sampleinfo(:,1));

% after sorting it should be the same
assert(isequal(dataAll.trialinfo(order1,:),  dataAppend.trialinfo(order2,:)));
assert(isequal(dataAll.sampleinfo(order1,:), dataAppend.sampleinfo(order2,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% timelock data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.keeptrials = 'yes';

timelockFIC = ft_timelockanalysis(cfg, dataFIC);
timelockFC  = ft_timelockanalysis(cfg, dataFC);
timelockAll = ft_timelockanalysis(cfg, dataAll);

cfg = [];
timelockAppend = ft_appendtimelock(cfg, timelockFIC, timelockFC);

%% compare the two versions

% these should not have sampleinfo as per specification in FT_DATATYPE_TIMELOCK
assert(~isfield(timelockAppend, 'sampleinfo'));
assert(~isfield(timelockAll,    'sampleinfo'));

% the order should be different
assert(~isequal(timelockAll.trialinfo,  timelockAppend.trialinfo));

% after sorting it should be the same
assert(isequal(timelockAll.trialinfo(order1,:),  timelockAppend.trialinfo(order2,:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% freq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.keeptrials = 'yes';
cfg.method = 'wavelet';
cfg.foilim = [5 15];
cfg.toi = -0.5:0.05:1.5;
cfg.pad = 'nextpow2';

freqFIC = ft_freqanalysis(cfg, dataFIC);
freqFC  = ft_freqanalysis(cfg, dataFC);
freqAll = ft_freqanalysis(cfg, dataAll);

cfg = [];
freqAppend = ft_appendfreq(cfg, freqFIC, freqFC);

%% compare the two versions

% these should not have sampleinfo as per specification in FT_DATATYPE_FREQ
assert(~isfield(freqAppend, 'sampleinfo'));
assert(~isfield(freqAll,    'sampleinfo'));

% the order should be different
assert(~isequal(freqAll.trialinfo,  freqAppend.trialinfo));

% after sorting it should be the same
assert(isequal(freqAll.trialinfo(order1,:),  freqAppend.trialinfo(order2,:)));
