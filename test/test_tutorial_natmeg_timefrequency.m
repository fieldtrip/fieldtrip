function test_tutorial_natmeg_timefrequency

% WALLTIME 00:30:00
% MEM 4gb

% this script executes the MATLAB content from
% http://www.fieldtriptoolbox.org/tutorial/natmeg/timefrequency
%
% it corresponds to the wiki version of 7 October 2014

clear all
close all

cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/natmeg'));

cfg = [];
cfg.dataset = 'oddball1_mc_downsampled.fif';
cfg.channel = 'MEG';

% define trials based on responses
cfg.trialdef.prestim       = 1.5;
cfg.trialdef.poststim      = 2.0;
cfg.trialdef.stim_triggers = [1 2];
cfg.trialdef.rsp_triggers  = [256 4096];
cfg.trialfun               = 'trialfun_oddball_responselocked';
cfg                        = ft_definetrial(cfg);

% preprocess MEG data
cfg.continuous             = 'yes';
cfg.demean                 = 'yes';
cfg.dftfilter              = 'yes';
cfg.dftfreq                = [50 100];

data_MEG_responselocked    = ft_preprocessing(cfg);

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.toi          = [-1 : 0.10 : 1.5];
cfg.foi          = 1:40;
cfg.t_ftimwin    = ones(size(cfg.foi)) * 0.5;

cfg.trials       = find(data_MEG_responselocked.trialinfo(:,1) == 256);
TFR_left_MEG     = ft_freqanalysis(cfg, data_MEG_responselocked);

cfg.trials       = find(data_MEG_responselocked.trialinfo(:,1) == 4096);
TFR_right_MEG    = ft_freqanalysis(cfg, data_MEG_responselocked);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.zlim         = [-2e-26 2e-26];
cfg.showlabels   = 'yes';
cfg.layout       = 'neuromag306mag.lay';
cfg.channel      = 'MEG*1';

figure;
ft_multiplotTFR(cfg, TFR_left_MEG);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.maskstyle    = 'saturation';
cfg.zlim         = [-1e-26 1e-26];
cfg.channel      = 'MEG1041';

figure;
ft_singleplotTFR(cfg, TFR_left_MEG);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.4 0.8];
cfg.zlim         = [-4e-27 4e-27];
cfg.ylim         = [15 25];
cfg.marker       = 'on';
cfg.layout       = 'neuromag306mag.lay';
cfg.channel      = 'MEG*1';

figure;
ft_topoplotTFR(cfg, TFR_left_MEG);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.4 0.8];
cfg.zlim         = [-4e-27 4e-27];
cfg.ylim         = [15 25];
cfg.marker       = 'on';
cfg.layout       = 'neuromag306mag.lay';
cfg.channel      = 'MEG*1';

figure;
ft_topoplotTFR(cfg, TFR_right_MEG);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)/(x1+x2)';

TFR_diff_MEG = ft_math(cfg, TFR_right_MEG, TFR_left_MEG);

cfg = [];
cfg.xlim         = [0.4 0.8];
cfg.zlim         = [-0.4 0.4];
cfg.ylim         = [15 25];
cfg.marker       = 'on';
cfg.layout       = 'neuromag306mag.lay';
cfg.channel      = 'MEG*1';

figure;
ft_topoplotTFR(cfg, TFR_diff_MEG);

cfg = [];
cfg.dataset = 'oddball1_mc_downsampled.fif';

% define trials based on responses
cfg.trialdef.prestim       = 1.5;
cfg.trialdef.poststim      = 2.0;
cfg.trialdef.stim_triggers = [1 2];
cfg.trialdef.rsp_triggers  = [256 4096];
cfg.trialfun               = 'trialfun_oddball_responselocked';
cfg                        = ft_definetrial(cfg);

% preprocess EEG data
cfg.channel                = 'EEG';
cfg.continuous             = 'yes';
cfg.demean                 = 'yes';
cfg.dftfilter              = 'yes';
cfg.dftfreq                = [50 100];

data_EEG_responselocked    = ft_preprocessing(cfg);

if false
  % skip the interactive section
  % select bad channels
  cfg = [];
  cfg.metric  = 'var';
  temp        = ft_rejectvisual(cfg, data_EEG_responselocked);
  % with this little trick we get the names of the selected channels
  badchannels = setdiff(data_EEG_responselocked.label,temp.label);
else
  badchannels = {
    'EEG001'
    'EEG002'
    'EEG003'
    'EEG004'
    'EEG005'
    'EEG006'
    'EEG007'
    'EEG008'
    'EEG015'
    'EEG016'
    'EEG017'
    'EEG097'
    'EEG098'
    'EEG099'
    'EEG100'
    'EEG101'
    'EEG102'
    'EEG111'
    'EEG112'};
end

% determine neighbours structure
cfg            = [];
cfg.method     = 'triangulation';
cfg.senstype   = 'EEG'; % Our data still contains information from the MEG channels, we want to make sure ft_prepare_neighbours does not get confused
neighbours_EEG = ft_prepare_neighbours(cfg, data_EEG_responselocked);

% plotting neighbours
cfg            = [];
cfg.neighbours = neighbours_EEG;
cfg.senstype   = 'EEG';
ft_neighbourplot(cfg, data_EEG_responselocked);

% fix channels
cfg = [];
cfg.method                    = 'spline';
cfg.neighbours                = neighbours_EEG;
cfg.badchannel                = badchannels;
cfg.senstype                  = 'EEG';
data_clean_EEG_responselocked = ft_channelrepair(cfg, data_EEG_responselocked);


cfg = [];
cfg.reref                  = 'yes';
cfg.refchannel             = 'all';
data_clean_EEG_responselocked = ft_preprocessing(cfg, data_clean_EEG_responselocked);

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.toi          = [-1 : 0.10 : 1.5];
cfg.foi          = 1:40;
cfg.t_ftimwin    = ones(size(cfg.foi)) * 0.5;

cfg.trials       = find(data_clean_EEG_responselocked .trialinfo(:,1) == 256);
TFR_left_EEG     = ft_freqanalysis(cfg, data_clean_EEG_responselocked );

cfg.trials       = find(data_clean_EEG_responselocked .trialinfo(:,1) == 4096);
TFR_right_EEG    = ft_freqanalysis(cfg, data_clean_EEG_responselocked );

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.5 1.0];
cfg.zlim         = [-4e-12 4e-12];
cfg.ylim         = [15 25];
cfg.marker       = 'on';
cfg.layout       = 'neuromag306eeg1005_natmeg.lay';

figure;
ft_topoplotTFR(cfg, TFR_left_EEG);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'relchange';
cfg.ylim         = [15 25];
cfg.xlim         = [0.5 1.0];
cfg.zlim         = [-1.2 1.2];
cfg.layout       = 'neuromag306eeg1005_natmeg.lay';

figure;
ft_topoplotTFR(cfg, TFR_left_EEG);

cfg = [];
cfg.parameter    = 'powspctrm';
cfg.operation    = '(x1-x2)/(x1+x2)';
TFR_diff_EEG = ft_math(cfg, TFR_right_EEG, TFR_left_EEG);

if false
  % skip the interactive section
  % if ft_math didn't work, then just do it by hand - its exactly the same:
  TFR_diff_EEG = TFR_right_EEG;
  TFR_diff_EEG.powspctrm = (TFR_right_EEG.powspctrm - TFR_left_EEG.powspctrm) ./ (TFR_right_EEG.powspctrm + TFR_left_EEG.powspctrm);
end

cfg = [];
cfg.xlim         = [0.4 0.8];
cfg.ylim         = [15 25];
cfg.zlim         = [-0.2 0.2];
cfg.marker       = 'on';
cfg.layout       = 'neuromag306eeg1005_natmeg.lay';

figure;
ft_topoplotTFR(cfg, TFR_diff_EEG);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.4 0.8];
cfg.ylim         = [15 25];
cfg.zlim         = [-1e-24 1e-24];
cfg.marker       = 'on';
cfg.layout       = 'neuromag306planar.lay';

figure;
ft_topoplotTFR(cfg, TFR_left_MEG);

TFR_left_MEG_comb  = ft_combineplanar([],TFR_left_MEG);
TFR_right_MEG_comb = ft_combineplanar([],TFR_right_MEG);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.4 0.8];
cfg.ylim         = [15 25];
cfg.zlim         = [-4e-24 4e-24];
cfg.marker       = 'on';
cfg.layout       = 'neuromag306cmb.lay';

figure;
ft_topoplotTFR(cfg, TFR_left_MEG_comb);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)/(x1+x2)';

TFR_diff_MEG_comb = ft_math(cfg, TFR_right_MEG_comb, TFR_left_MEG_comb);

cfg = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';
cfg.xlim         = [0.4 0.8];
cfg.ylim         = [15 25];
cfg.zlim         = [-0.3 0.3];
cfg.marker       = 'on';
cfg.layout       = 'neuromag306cmb.lay';

figure;
ft_topoplotTFR(cfg, TFR_diff_MEG_comb);
