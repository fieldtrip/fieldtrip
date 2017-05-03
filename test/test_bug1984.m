function test_bug1984

% WALLTIME 00:10:00
% MEM 3gb

dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.dataset                 = dataset;
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 2;
cfg.trialdef.eventvalue     = 3;                    % trigger value for fully incongruent (FIC)
cfg = ft_definetrial(cfg);
dataFIC = ft_preprocessing(cfg);

cfg = [];
cfg.dataset                 = dataset;
cfg.trialdef.eventtype      = 'backpanel trigger';
cfg.trialdef.prestim        = 1;
cfg.trialdef.poststim       = 1;                    % make these trials shorter
cfg.trialdef.eventvalue     = 9;                    % trigger value for fully congruent (FC)
cfg = ft_definetrial(cfg);
dataFC = ft_preprocessing(cfg);

cfg = [];
cfg.dataset = dataset;
cfg.trl = [
  dataFIC.cfg.trl
  dataFC.cfg.trl
  ];
[dum, order] = sort(cfg.trl(:,1));
cfg.trl = cfg.trl(order,:); % keep them in the original order of the experiment
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
%% freq data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.keeptrials = 'yes';
cfg.method = 'mtmfft';
cfg.foi = 1:10;
cfg.tapsmofrq = 2*(cfg.foi);

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
assert(~isequal(freqAll.cumsumcnt,  freqAppend.cumsumcnt));
assert(~isequal(freqAll.cumtapcnt,  freqAppend.cumtapcnt));

% after sorting it should be the same
assert(isequal(freqAll.trialinfo(order1,:),  freqAppend.trialinfo(order2,:)));
assert(isequal(freqAll.cumsumcnt(order1,:),  freqAppend.cumsumcnt(order2,:)));
assert(isequal(freqAll.cumtapcnt(order1,:),  freqAppend.cumtapcnt(order2,:)));


