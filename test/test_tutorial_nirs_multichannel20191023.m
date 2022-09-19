function test_tutorial_nirs_multichannel20191023

% WALLTIME 00:10:00
% MEM 4gb
% DEPENDENCY ft_datatype_sens ft_nirs_transform_ODs ft_nirs_prepare_ODtransformation

%%
% this reflects the "Preprocessing and averaging of multi-channel NIRS data" tutorial
% obtained from http://www.fieldtriptoolbox.org/tutorial/nirs_multichannel/ at 23 October 2019
% most of the comments have been removed to make this MATLAB script easier to read

cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/nirs_multichannel'));

%% ## Read data and downsample

cfg             = [];
cfg.dataset     = 'LR-01-2015-06-01-0002.oxy3';
data_raw        = ft_preprocessing(cfg);

cfg      = [];
cfg.opto = 'LR-01-2015-06-01-0002.oxy3';
ft_layoutplot(cfg);

%% ### Trigger channels

find(strcmp(data_raw.label,'ADC001'))
find(strcmp(data_raw.label,'ADC002'))

figure; hold on

plot(data_raw.time{1}, data_raw.trial{1}(97,:)*1.0, 'b-')
plot(data_raw.time{1}, data_raw.trial{1}(98,:)*1.1, 'r:')

event = ft_read_event('LR-01-2015-06-01-0002.oxy3');

data_raw.fsample

cfg                   = [];
cfg.resamplefs        = 10;
data_down             = ft_resampledata(cfg, data_raw);

cfg                = [];
cfg.preproc.demean = 'yes';
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.003   0.003 ];
cfg.channel        = 'Rx*'; % only show channels starting with Rx
ft_databrowser(cfg, data_down);

cfg                 = [];
cfg.hpfilter        = 'yes';
cfg.hpfreq          = 0.01;
data_flt            = ft_preprocessing(cfg,data_down);

cfg                = [];
cfg.preproc.demean = 'yes';
cfg.viewmode       = 'vertical';
cfg.continuous     = 'no';
cfg.ylim           = [ -0.003   0.003 ];
cfg.channel        = 'Rx*'; % only show channels starting with Rx
ft_databrowser(cfg, data_flt);

%% ## Epoch

event = ft_read_event('LR-01-2015-06-01-0002.oxy3');

adc001 = find(strcmp({event.type}, 'ADC001'));
adc002 = find(strcmp({event.type}, 'ADC002'));

% get the sample number in the original data
% note that we transpose them to get columns
smp001 = [event(adc001).sample]';
smp002 = [event(adc002).sample]';

factor = data_raw.fsample / data_down.fsample;

% get the sample number after downsampling
smp001 = round((smp001-1)/factor + 1);
smp002 = round((smp002-1)/factor + 1);

pre    =  round( 5*data_down.fsample);
post   =  round(20*data_down.fsample);
offset = -pre; % see ft_definetrial

trl001 = [smp001-pre smp001+post];
trl002 = [smp002-pre smp002+post];

% add the offset
trl001(:,3) = offset;
trl002(:,3) = offset;

trl001(:,4) = 1; % add a column with the condition number
trl002(:,4) = 2; % add a column with the condition number

% concatenate the two conditions and sort them
trl = sortrows([trl001; trl002]);

% remove trials that stretch beyond the end of the recording
sel = trl(:,2)<size(data_down.trial{1},2);
trl = trl(sel,:);

cfg     = [];
cfg.trl = trl;
data_epoch = ft_redefinetrial(cfg,data_down);

idx = find(data_epoch.trialinfo==2, 1, 'first')

cfg          = [];
cfg.channel  = 'Rx*';
cfg.trials   = 8;
cfg.baseline = 'yes';
ft_singleplotER(cfg, data_epoch);

%% ## Remove bad channels

cfg                 = [];
data_sci            = ft_nirs_scalpcouplingindex(cfg, data_epoch);

%% ## Transform optical densities to oxy- and deoxy-hemoglobin concentration changes

cfg                 = [];
cfg.target          = {'O2Hb', 'HHb'};
cfg.channel         = 'nirs'; % e.g. one channel incl. wildcards, you can also use 'all' to select all nirs channels
data_conc           = ft_nirs_transform_ODs(cfg, data_sci);

%% ## Separate functional from systemic responses
%
%% ### Low-pass filtering

cfg                   = [];
cfg.lpfilter          = 'yes';
cfg.lpfreq            = 0.8;
data_lpf              = ft_preprocessing(cfg, data_conc);

%% ## Plot results

cfg               = [];
cfg.trials        = find(data_lpf.trialinfo(:,1) == 1);
timelockSTD       = ft_timelockanalysis(cfg, data_lpf);

cfg                 = [];
cfg.baseline        = [-5 0];
timelockSTD         = ft_timelockbaseline(cfg, timelockSTD);

cfg           = [];
cfg.trials    = find(data_lpf.trialinfo(:,1) == 2);
timelockDEV   = ft_timelockanalysis(cfg, data_lpf);

cfg           = [];
cfg.baseline  = [-5 0];
timelockDEV   = ft_timelockbaseline(cfg, timelockDEV);

%%
% The lines below are not exactly according to the 20191023 version, but have been updated to make it work.
% The tutorial on the website has also been updated, see https://github.com/fieldtrip/fieldtrip/pull/1245#issuecomment-546270124

load('nirs_48ch_layout.mat')
figure; ft_plot_layout(lay) % note that O2Hb and HHb channels fall on top of each other

cfg                   = [];
cfg.showlabels        = 'yes';
cfg.layout            = lay;
cfg.interactive       = 'yes';
cfg.graphcolor        = 'rb';
cfg.colorgroups       = ones(numel(timelockDEV.label),1);
cfg.colorgroups(contains(timelockDEV.label, 'O2Hb')) = 1; % these will be red
cfg.colorgroups(contains(timelockDEV.label, 'HHb'))  = 2; % these will be blue
ft_multiplotER(cfg, timelockDEV);

cfg          = [];
cfg.layout   = lay;
cfg.marker   = 'labels';
cfg.xlim     = [5 7];
cfg.zlim     = [-0.2 0.2];
cfg.channel  = '* [O2Hb]';
figure; subplot(1,2,1);
ft_topoplotER(cfg, timelockDEV);
title('[O2Hb]');

cfg.channel  = '* [HHb]';
subplot(1,2,2);
ft_topoplotER(cfg, timelockDEV);
title('[HHb]');
