function test_ft_redefinetrial

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY

%% use 10 trials from the ctf151 data structure
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'));

data.trialinfo = (1:10)';

cfg        = [];
cfg.trials = 1:2:10;
data2      = ft_redefinetrial(cfg, data);
assert(all(data2.trialinfo==[1 3 5 7 9]'));

cfg          = [];
cfg.trl      = data.sampleinfo;
cfg.trl(:,1) = cfg.trl(:,1)+50;
cfg.trl(:,3) = -50;
data3        = ft_redefinetrial(cfg, data);
assert(all(data3.sampleinfo(:,1)==[51:300:2751]'));

cfg        = [];
cfg.toilim = [0.2 0.8];
data4      = ft_redefinetrial(cfg, data);
assert(all(data4.sampleinfo(:,1)==[61:300:2761]') && all(data4.sampleinfo(:,2)==[241:300:2941]'));

cfg        = [];
cfg.length = 0.2;
data5      = ft_redefinetrial(cfg, data);
assert(numel(data5.trial)==50 && all(data5.trialinfo(1:5)==1));

% this should give an error
% cfg        = [];
% cfg.trl    = [1 600 0];
% data6      = ft_redefinetrial(cfg, data);
% assert(~isfield(data6, 'trialinfo'));

% here the trialinfo should be added
data = rmfield(data, 'trialinfo');
cfg          = [];
cfg.trl      = data.sampleinfo;
cfg.trl(:,1) = cfg.trl(:,1)+50;
cfg.trl(:,3) = -50;
cfg.trl(:,4) = 1:10';
data7        = ft_redefinetrial(cfg, data);
assert(all(data7.trialinfo(:,1)==(1:10)'));

% here the user-specified trialinfo should be added
data.trialinfo = (1:10)';
cfg          = [];
cfg.trl      = data.sampleinfo;
cfg.trl(:,1) = cfg.trl(:,1)+50;
cfg.trl(:,3) = -50;
cfg.trl(:,4) = (11:20)';
data8        = ft_redefinetrial(cfg, data);
assert(all(data8.trialinfo(:,1)==(11:20)'));

% here the old trialinfo should be added
data.trialinfo = (1:10)';
cfg          = [];
cfg.trl      = data.sampleinfo;
cfg.trl(:,1) = cfg.trl(:,1)+50;
cfg.trl(:,3) = -50;
data8        = ft_redefinetrial(cfg, data);
assert(all(data8.trialinfo(:,1)==(1:10)'));

% we can specify separate toilim per trial
cfg        = [];
cfg.toilim = repmat([0.2 0.8], numel(data.trial), 1);
data9      = ft_redefinetrial(cfg, data);
assert(all(data9.sampleinfo(:,1)==[61:300:2761]') && all(data9.sampleinfo(:,2)==[241:300:2941]'));

%% construct a continuous data structure

data_orig = [];
data_orig.label = {'1'};
data_orig.time{1} = ((1:10000)-1)./1000; % 10 seconds
data_orig.trial{1} = randn(1, 10000);
data_orig.sampleinfo = [1 100000];
data_orig.trialinfo = [1];

cfg = [];
cfg.length = 1;
data_segmented = ft_redefinetrial(cfg, data_orig);
assert(length(data_segmented.trial)==10);

cfg = [];
cfg.continuous = 'yes';
data_continuous = ft_redefinetrial(cfg, data_segmented);
assert(length(data_continuous.trial)==1);

cfg = [];
cfg.trials = setdiff(1:10, 5); % remove one trial
data_segmented = ft_selectdata(cfg, data_segmented);

cfg = [];
cfg.continuous = 'yes';
data_continuous = ft_redefinetrial(cfg, data_segmented);
assert(length(data_continuous.trial)==2);
