function test_tutorial_sensor_analysis(datadir)

% MEM 2500mb
% WALLTIME 00:25:00

% TEST test_tutorial_sensor_overview
% TEST ft_redefinetrial ft_freqanalysis ft_timelockanalysis ft_appenddata ft_prepare_neighbours ft_megplanar ft_combineplanar ft_multiplotER ft_multiplotTFR ft_connectivityanalysis

global ft_default;
ft_default.feedback = 'no';
ft_default.trackconfig = 'no'; % don't convert the cfg into a config object, as that fails in r8540 due to a subsref error

if nargin==0
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/sensor_analysis');
end

load(fullfile(datadir, 'subjectK.mat'));

%% show a trial

plot(data_left.time{1}, data_left.trial{1}(130,:));

for k = 1:10
  plot(data_left.time{k}, data_left.trial{k}(130,:)+k*1.5e-12);
  hold on;
end
plot([0 0], [0 1], 'k');
ylim([0 11*1.5e-12]);
set(gca, 'ytick', (1:10).*1.5e-12);
set(gca, 'yticklabel', 1:10);
ylabel('trial number');
xlabel('time (s)');

%% event-related field

% combine the data
data = ft_appenddata([], data_left, data_right);

% timelock analysis
cfg                 = [];
cfg.channel         = 'MEG';
cfg.vartrllength    = 1;
tl                  = ft_timelockanalysis(cfg, data);

% plot
cfg                 = [];
cfg.showlabels      = 'yes';
cfg.showoutline     = 'yes';
cfg.layout          = 'CTF151.lay';
ft_multiplotER(cfg, tl);

% repeat with planar gradient

cfg                 = [];
cfg.method          = 'template';
cfg.template        = 'CTF151_neighb.mat';
cfg.feedback        = 'yes';
neighbours          = ft_prepare_neighbours(cfg, data);

cfg                 = [];
cfg.method          = 'sincos';
cfg.neighbours      = neighbours;
data_planar         = ft_megplanar(cfg, data);

% timelock analysis
cfg                 = [];
cfg.channel         = 'MEG';
cfg.vartrllength    = 1;
tl_planar           = ft_timelockanalysis(cfg, data_planar);

cfg                 = [];
tl_plancmb          = ft_combineplanar(cfg, tl_planar);

% plot
cfg                 = [];
cfg.showlabels      = 'yes';
cfg.showoutline     = 'yes';
cfg.layout          = 'CTF151.lay';
ft_multiplotER(cfg, tl_plancmb);

%% time-frequency analysis

% subselect data to speed things up
cfg                 = [];
cfg.toilim          = [-0.8 1];
cfg.minlength       = 'maxperlen'; % this ensures all trials are equal length
data_small          = ft_redefinetrial(cfg, data_planar);

cfg                 = [];
cfg.method          = 'mtmconvol';
cfg.taper           = 'hanning';
cfg.channel         = 'MEG';

% set the frequencies of interest
cfg.foi             = 20:5:100;

% set the timepoints of interest: from -0.8 to 1.1 in steps of 100ms
cfg.toi             = -0.8:0.1:1;

% set the time window for TFR analysis: constant length of 200ms
cfg.t_ftimwin       = 0.2 * ones(length(cfg.foi), 1);

% average over trials
cfg.keeptrials      = 'no';

% pad trials to integer number of seconds, this speeds up the analysis
% and results in a neatly spaced frequency axis
cfg.pad             = 2;
freq                = ft_freqanalysis(cfg, data_small);

cfg                 = [];
freq                = ft_combineplanar(cfg, freq);

% plot
cfg                 = [];
cfg.interactive     = 'yes';
cfg.showoutline     = 'yes';
cfg.layout          = 'CTF151.lay';
cfg.baseline        = [-0.8 0];
cfg.baselinetype    = 'relchange';
cfg.zlim            = 'maxabs';
ft_multiplotTFR(cfg, freq);

%% sensor-level coherence

% show an EMG trial
figure
subplot(2,1,1);
plot(data_left.time{12},data_left.trial{12}(77,:));
axis tight;
legend(data_left.label(77));

subplot(2,1,2);
plot(data_left.time{12},data_left.trial{12}(152:153,:));
axis tight;
legend(data_left.label(152:153));

% subselect stimulation period
cfg                 = [];
cfg.toilim          = [-1 -0.0025];
cfg.minlength       = 'maxperlen'; % this ensures all trials are equal length
data_stim           = ft_redefinetrial(cfg, data);

% estimate power and csd
cfg                 = [];
cfg.output          = 'powandcsd';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.foilim          = [5 100];
cfg.tapsmofrq       = 5;
cfg.keeptrials      = 'yes';
cfg.channel         = {'MEG' 'EMGlft' 'EMGrgt'};
cfg.channelcmb      = {'MEG' 'EMGlft'; 'MEG' 'EMGrgt'};
freq_csd            = ft_freqanalysis(cfg, data_stim);

% connectivity analysis
cfg                 = [];
cfg.method          = 'coh';
cfg.channelcmb      = {'MEG' 'EMG'};
conn                = ft_connectivityanalysis(cfg, freq_csd);

% plot
cfg                 = [];
cfg.parameter       = 'cohspctrm';
cfg.xlim            = [5 80];
cfg.refchannel      = 'EMGrgt';
cfg.layout          = 'CTF151.lay';
cfg.showlabels      = 'no';
figure;
ft_multiplotER(cfg, conn);

end
