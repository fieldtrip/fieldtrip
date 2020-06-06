function test_issue1410

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'), 'data');
grad = data.grad;
clear data

%%

vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.unit = 'cm';
vol.type = 'singlesphere';

%%

cfg = [];
cfg.dip.pos = [0 0 9];
cfg.dip.mom = [1 0 0];
cfg.dip.unit = 'cm';
cfg.headmodel = vol;
cfg.grad = grad;
data = ft_dipolesimulation(cfg);

cfg = [];
timelock = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'fourier';
freq = ft_freqanalysis(cfg, data);

%%

tfdata = timelock;
method = 'lcmv';

cfg = [];
cfg.xgrid = -8:2:8;
cfg.ygrid = -8:2:8;
cfg.zgrid = (-8:2:8) + 4;
cfg.unit = 'cm';
cfg.headmodel = vol;
cfg.method = method;
cfg.keepleadfield = 'yes';
cfg.(method).keepfilter = 'yes';
cfg.(method).lambda = '10%';
cfg.channel = 'MEG';
source1 = ft_sourceanalysis(cfg, tfdata);

%%

% in the next code the labels are missing for the precomputed filters
% this should assume that they were computed with the same channel selection

cfg = [];
cfg.method = method;
cfg.sourcemodel.pos = source1.pos;
cfg.sourcemodel.filter = source1.avg.filter;
cfg.channel = 'MEG';
source2 = ft_sourceanalysis(cfg, tfdata);
