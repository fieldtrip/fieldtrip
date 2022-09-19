function test_pull1138

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourceplot ft_plot_cloud

load('/home/common/matlab/fieldtrip/data/ftp/tutorial/human_ecog/SubjectUCI29/SubjectUCI29_freq.mat');
cortex = load('/home/common/matlab/fieldtrip/data/ftp/tutorial/human_ecog/SubjectUCI29/SubjectUCI29_hull_lh.mat');

%%

cfg = [];
cfg.baseline = [-.3 -.1];
cfg.baselinetype = 'relchange';
freq_blc = ft_freqbaseline(cfg, freq);

cfg = [];
cfg.frequency = [70 150];
cfg.avgoverfreq = 'yes';
cfg.latency = [0 0.8];
cfg.avgovertime = 'yes';
freq_sel = ft_selectdata(cfg, freq_blc);

%% default
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud'; % <-- uses ft_plot_cloud
ft_sourceplot(cfg, freq_sel, cortex.mesh);
view([120 40]);
lighting gouraud;
camlight;

%% cloudtype cloud
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud';
cfg.cloudtype = 'cloud';  % <-- uses ft_plot_cloud
ft_sourceplot(cfg, freq_sel, cortex.mesh);
view([120 40]);
lighting gouraud;
camlight;

%% cloudtype surf
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud';
cfg.cloudtype = 'surf';  % <-- uses ft_plot_cloud
ft_sourceplot(cfg, freq_sel, cortex.mesh);
view([120 40]);
lighting gouraud;
camlight;

%% cloudtype point
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud';
cfg.cloudtype = 'point';  % <-- uses plot3
ft_sourceplot(cfg, freq_sel, cortex.mesh);
view([120 40]);
lighting gouraud;
camlight;

