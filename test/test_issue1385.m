function test_issue1385

% MEM 4gb
% WALLTIME 00:20:00
% DEPENDENCY ft_sourceplot

datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/human_ecog');
subj = 'SubjectUCI29';
subj_dir = fullfile(datadir, subj);
cd(subj_dir)

load('SubjectUCI29_freq.mat')
pial_lh = ft_read_headshape('freesurfer/surf/lh.pial');

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

cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'surface';
cfg.interpmethod = 'sphere_weighteddistance';
cfg.sphereradius = 8;
cfg.camlight = 'no';
ft_sourceplot(cfg, freq_sel, pial_lh);