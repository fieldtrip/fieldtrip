function test_tutorial_natmeg_statistics

% WALLTIME 00:20:00
% MEM 4gb

% this script executes the MATLAB content from
% http://www.fieldtriptoolbox.org/tutorial/natmeg/statistics
%
% it corresponds to the wiki version of 7 October 2014

clear all
close all

cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/natmeg'));

cfg = [];
cfg.dataset = 'oddball1_mc_downsampled.fif';
% define trials based on response
cfg.trialdef.prestim       = 1.0;
cfg.trialdef.poststim      = 2.0;
cfg.trialdef.stim_triggers = [1 2];
cfg.trialdef.rsp_triggers  = [256 4096];
cfg.trialfun               = 'trialfun_oddball_responselocked';
cfg                        = ft_definetrial(cfg);
% preprocess MEG data
cfg.channel = 'MEG*1';
cfg.continuous             = 'yes';
cfg.demean                 = 'yes';
cfg.dftfilter              = 'yes';
cfg.dftfreq                = [50 100];
data_responselocked        = ft_preprocessing(cfg);

cfg           = [];
cfg.output    = 'pow';
cfg.method    = 'mtmconvol';
cfg.taper     = 'hanning';
cfg.toi       = 0.0 : 0.1 : 1.0;
cfg.foi       = 12:24;
cfg.t_ftimwin = ones(size(cfg.foi)) * 0.5;

cfg.trials    = find(data_responselocked.trialinfo(:,1) == 256);
TFR_left      = ft_freqanalysis(cfg, data_responselocked);

cfg.trials    = find(data_responselocked.trialinfo(:,1) == 4096);
TFR_right     = ft_freqanalysis(cfg, data_responselocked);
cfg = [];
cfg.baseline     = [0 0];
cfg.baselinetype = 'absolute';
cfg.layout       = 'neuromag306mag.lay';

figure;
ft_multiplotTFR(cfg, TFR_left);

figure;
ft_multiplotTFR(cfg, TFR_right);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = '(x1-x2)/(x1+x2)';
TFR_diff      = ft_math(cfg, TFR_right, TFR_left);

cfg = [];
cfg.marker  = 'on';
cfg.layout  = 'neuromag306mag.lay';
cfg.channel = 'MEG*1';
figure; ft_multiplotTFR(cfg, TFR_diff);

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.toi          = 0.0 : 0.1 : 1.0;
cfg.foi          = 15:25;
cfg.t_ftimwin    = ones(size(cfg.foi)) * 0.5;
cfg.keeptrials   = 'yes';  % keep the TFR on individual trials
TFR_all          = ft_freqanalysis(cfg, data_responselocked);

cfg           = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
TFR_logpow    = ft_math(cfg, TFR_all);


cfg           = [];
cfg.channel   = 'MEG*1';
cfg.method    = 'triangulation';
cfg.grad      = TFR_all.grad;
cfg.feedback  = 'yes';
neighbours    = ft_prepare_neighbours(cfg);

cfg = [];
cfg.channel   = 'MEG*1';
cfg.statistic = 'indepsamplesT';
cfg.ivar      = 1;
cfg.design    = zeros(1, size(TFR_all.trialinfo,1));

cfg.design(TFR_all.trialinfo(:,1)== 256) = 1;
cfg.design(TFR_all.trialinfo(:,1)==4096) = 2;

cfg.method    = 'analytic';
cfg.correctm  = 'no';
TFR_stat1     = ft_freqstatistics(cfg, TFR_logpow);

cfg.method    = 'analytic';
cfg.correctm  = 'bonferoni';
TFR_stat2     = ft_freqstatistics(cfg, TFR_logpow);
cfg.method    = 'analytic';
cfg.correctm  = 'fdr';
TFR_stat3     = ft_freqstatistics(cfg, TFR_logpow);
cfg.method            = 'montecarlo';
cfg.correctm          = 'cluster';
cfg.numrandomization  = 500; % 1000 is recommended, but takes longer
cfg.neighbours        = neighbours;
TFR_stat4     = ft_freqstatistics(cfg, TFR_logpow);

cfg               = [];
cfg.marker        = 'on';
cfg.layout        = 'neuromag306mag.lay';
cfg.channel       = 'MEG*1';
cfg.parameter     = 'stat';  % plot the t-value
cfg.maskparameter = 'mask';  % use the thresholded probability to mask the data
cfg.maskstyle     = 'saturation';
figure; ft_multiplotTFR(cfg, TFR_stat1);
figure; ft_multiplotTFR(cfg, TFR_stat2);
figure; ft_multiplotTFR(cfg, TFR_stat3);
figure; ft_multiplotTFR(cfg, TFR_stat4);

cfg = [];
cfg.dataset = 'oddball1_mc_downsampled.fif';

% define trials based on stimulus
cfg.trialdef.prestim  = 0.3;
cfg.trialdef.poststim = 0.7;
cfg.trialdef.stim_triggers = [1 2];
cfg.trialdef.rsp_triggers  = [256 4096];
cfg.trialfun          = 'trialfun_oddball_stimlocked';
cfg                   = ft_definetrial(cfg);

% preprocess MEG data
cfg.channel           = 'MEG*1';
cfg.continuous        = 'yes';
cfg.demean            = 'yes';
cfg.baselinewindow    = [-inf 0];
cfg.dftfilter         = 'yes';
cfg.dftfreq           = [50 100];
data_stimlocked = ft_preprocessing(cfg);

cfg                 = [];
cfg.trials          = 1:100;
data_stimlocked     = ft_selectdata(cfg, data_stimlocked);

cfg         = [];
cfg.trials  = find(data_stimlocked.trialinfo(:,1) == 1);
ERF_std     = ft_timelockanalysis(cfg, data_stimlocked);

cfg.trials  = find(data_stimlocked.trialinfo(:,1) == 2);
ERF_dev     = ft_timelockanalysis(cfg, data_stimlocked);

cfg             = [];
cfg.latency     = [0.08 0.11];
cfg.avgovertime = 'yes';
ERF_peak        =  ft_selectdata(cfg, ERF_std);

cfg        = [];
cfg.layout = 'neuromag306mag.lay';
figure; ft_multiplotER(cfg, ERF_std, ERF_dev);

cfg             = [];
cfg.keeptrials  = 'yes';
ERF_all = ft_timelockanalysis(cfg, data_stimlocked);

cfg           = [];
cfg.statistic = 'indepsamplesT';
cfg.design    = zeros(1, size(ERF_all.trialinfo,1));
cfg.ivar      = 1;
cfg.design(ERF_all.trialinfo(:,1)==1) = 1;
cfg.design(ERF_all.trialinfo(:,1)==2) = 2;
cfg.method    = 'analytic';
cfg.correctm  = 'no';
ERF_stat1     = ft_timelockstatistics(cfg, ERF_all);

cfg.method    = 'analytic';
cfg.correctm  = 'bonferoni';
ERF_stat2     = ft_timelockstatistics(cfg, ERF_all);

cfg.method    = 'analytic';
cfg.correctm  = 'fdr';
ERF_stat3     = ft_timelockstatistics(cfg, ERF_all);

cfg.method            = 'montecarlo';
cfg.correctm          = 'cluster';
cfg.numrandomization  = 500; % 1000 is recommended, but that takes longer
cfg.neighbours        = neighbours;
ERF_stat4     = ft_timelockstatistics(cfg, ERF_all);

cfg = [];
cfg.layout        = 'neuromag306mag.lay';
cfg.maskparameter = 'mask';
cfg.maskstyle     = 'box';
ERF_std.mask = ERF_stat1.mask;  % copy the significance mask into the ERF
figure; ft_multiplotER(cfg, ERF_std, ERF_dev);
title('no correction');

ERF_std.mask = ERF_stat2.mask;  % copy the significance mask into the ERF
figure; ft_multiplotER(cfg, ERF_std, ERF_dev);
title('bonferroni');

ERF_std.mask = ERF_stat3.mask;  % copy the significance mask into the ERF
figure; ft_multiplotER(cfg, ERF_std, ERF_dev);
title('fdr');

ERF_std.mask = ERF_stat4.mask;  % copy the significance mask into the ERF
figure; ft_multiplotER(cfg, ERF_std, ERF_dev);
title('cluster');
