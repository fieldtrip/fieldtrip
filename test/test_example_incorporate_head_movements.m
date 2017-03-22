function test_example_incorporate_head_movements

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_definetrial ft_preprocessing ft_timelockanalysis ft_regressconfound

global ft_default;
ft_default.feedback = 'no';

dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/example/regressconfound/TacStimRegressConfound.ds');

cfg                         = [];
cfg.dataset                 = dataset;
cfg.trialdef.eventtype      = 'UPPT001';
cfg.trialdef.eventvalue     = 4;
cfg.trialdef.prestim        = 0.2;
cfg.trialdef.poststim       = 0.3;
cfg.continuous              = 'yes';
cfg = ft_definetrial(cfg);

cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
  'HLC0021','HLC0022','HLC0023', ...
  'HLC0031','HLC0032','HLC0033'};

headpos = ft_preprocessing(cfg);

% calculate the mean coil position per trial
ntrials = length(headpos.sampleinfo);
for t = 1:ntrials
  coil1(:,t) = [mean(headpos.trial{1,t}(1,:)); mean(headpos.trial{1,t}(2,:)); mean(headpos.trial{1,t}(3,:))];
  coil2(:,t) = [mean(headpos.trial{1,t}(4,:)); mean(headpos.trial{1,t}(5,:)); mean(headpos.trial{1,t}(6,:))];
  coil3(:,t) = [mean(headpos.trial{1,t}(7,:)); mean(headpos.trial{1,t}(8,:)); mean(headpos.trial{1,t}(9,:))];
end

% calculate the headposition and orientation per trial (for function see bottom page)
cc = circumcenter(coil1, coil2, coil3);

cc_rel = [cc - repmat(cc(:,1),1,size(cc,2))]';

% plot translations
figure();
plot(cc_rel(:,1:3)*1000) % in mm

% plot rotations
figure();
plot(cc_rel(:,4:6))

maxposchange = max(abs(cc_rel(:,1:3)*1000)) % in mm;

% define trials
cfg                         = [];
cfg.dataset                 = dataset;
cfg.trialdef.eventtype      = 'UPPT001';
cfg.trialdef.eventvalue     = 4;
cfg.trialdef.prestim        = 0.2;
cfg.trialdef.poststim       = 0.3;
cfg.continuous              = 'yes';
cfg = ft_definetrial(cfg);

% preprocess the MEG data
cfg.channel                 = {'MEG'};
cfg.demean                  = 'yes';
cfg.baselinewindow          = [-0.2 0];
cfg.dftfilter               = 'yes'; % notch filter to filter out 50Hz
data = ft_preprocessing(cfg);

% timelock analysis
cfg                         = [];
cfg.keeptrials              = 'yes';
timelock = ft_timelockanalysis(cfg, data);

% define trials
cfg                         = [];
cfg.dataset                 = dataset;
cfg.trialdef.eventtype      = 'UPPT001';
cfg.trialdef.eventvalue     = 4;
cfg.trialdef.prestim        = 0.2;
cfg.trialdef.poststim       = 0.3;
cfg.continuous              = 'yes';
cfg = ft_definetrial(cfg);

% preprocess the headposition data
cfg.channel                 = {'HLC0011','HLC0012','HLC0013', ...
  'HLC0021','HLC0022','HLC0023', ...
  'HLC0031','HLC0032','HLC0033'};
headpos = ft_preprocessing(cfg);

% calculate the mean coil position per trial
ntrials = length(headpos.sampleinfo);
for t = 1:ntrials
  coil1(:,t) = [mean(headpos.trial{1,t}(1,:)); mean(headpos.trial{1,t}(2,:)); mean(headpos.trial{1,t}(3,:))];
  coil2(:,t) = [mean(headpos.trial{1,t}(4,:)); mean(headpos.trial{1,t}(5,:)); mean(headpos.trial{1,t}(6,:))];
  coil3(:,t) = [mean(headpos.trial{1,t}(7,:)); mean(headpos.trial{1,t}(8,:)); mean(headpos.trial{1,t}(9,:))];
end

% calculate the headposition and orientation per trial
cc = circumcenter(coil1, coil2, coil3);

% demean to obtain translations and rotations from the average position and orientation
cc_dem = [cc - repmat(mean(cc,2),1,size(cc,2))]';

% add head movements, but also the squares, cubes and all their derivatives to
% the regressorlist. also add the constant (at the end; column 37)
confound = [cc_dem cc_dem.^2 cc_dem.^3 ...
  gradient(cc_dem')' gradient(cc_dem.^2')' gradient(cc_dem.^3')' ...
  ones(size(cc_dem,1),1)];

% regress out headposition confounds
cfg                         = [];
cfg.confound                = confound;
cfg.reject                  = [1:36]; % keeping the constant (nr 37)
regr = ft_regressconfound(cfg, timelock);


