function test_issue1521

% MEM 4gb
% WALLTIME 00:10:00
% DEPENDENCY ft_combineplanar svdfft ft_dipolesimulation

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1521.mat'));

%%

% cfg.method = 'sum' only works for frequency data with powspctrm
% cfg.method = 'abssvd' is not supported for frequency data
% cfg.method = 'complex' is not supported for frequency data

cfg = [];
cfg.method = 'svd';
freqplanar_svd = ft_combineplanar(cfg, freqplanar);

cfg = [];
powplanar = ft_freqdescriptives(cfg, freqplanar);

cfg = [];
cfg.method = 'sum';
powplanar_sum = ft_combineplanar(cfg, powplanar);

%%

cfg = [];
cfg.headmodel.type = 'singlesphere';
cfg.headmodel.o = [0 0 4];
cfg.headmodel.r = 12;
cfg.headmodel.cond = 1;
cfg.headmodel.unit = 'cm';
cfg.channel = 'MEGGRAD';
cfg.grad = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds'), 'senstype', 'meg');
raw = ft_dipolesimulation(cfg);

cfg = [];
cfg.method = 'template';
nb = ft_prepare_neighbours(cfg, raw);

cfg = [];
cfg.method = 'sincos';
cfg.neighbours = nb;
rawplanar = ft_megplanar(cfg, raw);

cfg = [];
timelockplanar = ft_timelockanalysis(cfg, rawplanar);

%%

cfg = [];
cfg.method = 'svd';
timelockplanar_svd = ft_combineplanar(cfg, timelockplanar);

cfg = [];
cfg.method = 'abssvd';
timelockplanar_abssvd = ft_combineplanar(cfg, timelockplanar);

cfg = [];
cfg.method = 'sum';
timelockplanar_sum = ft_combineplanar(cfg, timelockplanar);

cfg = [];
cfg.method = 'complex';
timelockplanar_complex = ft_combineplanar(cfg, timelockplanar);
