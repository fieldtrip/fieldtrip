function test_tutorial_nirs_singlechannel20191023

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_datatype_sens ft_nirs_transform_ODs ft_nirs_prepare_ODtransformation

%%
% this reflects the "Preprocessing and averaging of single-channel NIRS data" tutorial
% obtained from http://www.fieldtriptoolbox.org/tutorial/nirs_singlechannel/ at 23 October 2019
% most of the comments have been removed to make this MATLAB script easier to read

cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/nirs_singlechannel'));


%% ## Read & trim data

cfg = [];
cfg.dataset = 'motor_cortex.oxy3';

[data] = ft_preprocessing(cfg);

cfg = [];
cfg.ylim = 'maxmin';
ft_databrowser(cfg, data);

cfg = [];
cfg.ylim = 'maxmin';
cfg.channel = {'Rx4b-Tx5 [860nm]', 'Rx4b-Tx5 [764nm]'};
ft_databrowser(cfg, data);

%% ## Remove artifacts

cfg = [];
cfg.channel = {'Rx4b-Tx5 [860nm]', 'Rx4b-Tx5 [764nm]'};
cfg.artfctdef.zvalue.channel = {'Rx4b-Tx5 [860nm]', 'Rx4b-Tx5 [764nm]'};
cfg.artfctdef.zvalue.cutoff = 3.5;
[cfg, artifact] = ft_artifact_zvalue(cfg, data);

%% ## Transform to changes in oxyHB/deoxyHB

cfg = [];
cfg.dpf = 5.9;
cfg.channel = {'Rx4b-Tx5 [860nm]', 'Rx4b-Tx5 [764nm]'};
data_conc = ft_nirs_transform_ODs(cfg, data);

%% ## Separate functional from systemic responses

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [0.01 0.1];
data_filtered = ft_preprocessing(cfg, data_conc);

%% ## Define epochs of interest

cfg = [];
cfg.dataset = 'motor_cortex.oxy3';
cfg.channel = {'Rx4b-Tx5 [860nm]', 'Rx4b-Tx5 [764nm]'};
cfg.trialdef = [];
cfg.trialdef.eventtype = '?';

ft_definetrial(cfg);

cfg.trialdef.eventtype  = 'event';
cfg.trialdef.eventvalue = 'A';
cfg.trialdef.prestim    = 10;
cfg.trialdef.poststim   = 35;
cfg = ft_definetrial(cfg);
data_epoch = ft_redefinetrial(cfg,data_filtered);

%
cfg = [];
cfg.ylim = [-1 1];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, data_epoch);

%% ## Timelockanalysis

cfg = [];
data_timelock = ft_timelockanalysis(cfg, data_epoch);

time = data_timelock.time;
O2Hb = data_timelock.avg(1,:);
HHb  = data_timelock.avg(2,:);

figure;
plot(time,O2Hb,'r'); hold on;
plot(time,HHb,'b');
legend('O2Hb','HHb'); ylabel('\DeltaHb (\muM)'); xlabel('time (s)');

