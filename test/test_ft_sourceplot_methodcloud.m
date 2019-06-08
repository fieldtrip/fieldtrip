function test_ft_sourceplot_methodcloud

% MEM 1500mb
% WALLTIME 00:10:00

% FIXME: please change location to Donders location of human ecog tutorial dataset
load('/Users/arjsto/Documents/Talks/FieldTrip/davis2019/data/human_ecog/SubjectUCI29_freq.mat');
cortex = load('/Users/arjsto/Documents/Talks/FieldTrip/davis2019/data/human_ecog/SubjectUCI29_hull_lh.mat');

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

% default
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud';
ft_sourceplot(cfg, freq_sel, cortex.mesh);
view([120 40]);
lighting gouraud;
camlight;

% cloudtype cloud
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud';
cfg.cloudtype = 'cloud';  % <--
ft_sourceplot(cfg, freq_sel, cortex.mesh);
view([120 40]);
lighting gouraud;
camlight;

% cloudtype point
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud';
cfg.cloudtype = 'point';  % <--
ft_sourceplot(cfg, freq_sel, cortex.mesh);
view([120 40]);
lighting gouraud;
camlight;

% cloudtype surf
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud';
cfg.cloudtype = 'surf';  % <--
ft_sourceplot(cfg, freq_sel, cortex.mesh);
view([120 40]);
lighting gouraud;
camlight;
