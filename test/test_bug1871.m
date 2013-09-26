function test_bug1871

% WALLTIME 00:03:33

% TEST test_bug1871
% TEST ft_struct2single 

cd(dccnfilename('/home/common/matlab/fieldtrip/data/test'));
load avgFIC.mat

avgFIC.avg = ft_struct2single(avgFIC.avg);

cfg = [];
cfg.xlim = [0.3 0.5];
cfg.layout = 'CTF151.lay';
cfg.zlim = 'maxmin';
cfg.colorbar = 'yes';
ft_topoplotER(cfg,avgFIC);
ft_multiplotER(cfg, avgFIC);

load bug1871

freq.powspctrm = ft_struct2single(freq.powspctrm);

cfg = [];
cfg.channel     = 'MEGGRAD';
cfg.interactive = 'yes';
cfg.showlabels  = 'no';
cfg.zlim        = 'maxabs';
cfg.layout      = 'neuromag306cmb.lay';
ft_topoplotTFR(cfg, freq);
ft_multiplotTFR(cfg, freq);

