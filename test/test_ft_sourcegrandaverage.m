function test_ft_sourcegrandaverage

% WALLTIME 00:10:00
% MEM 1gb

% TEST test_bug2185
% TEST ft_sourcegrandaverage ft_selectdata ft_selectdata_new ft_datatype_source

global ft_default
ft_default = [];

%% format 1

source = [];
source.dim = [10 11 12];
source.transform = eye(4);
source.avg.pow = rand(10*11*12,1);
source.inside = 1:660;
source.outside = 661:1320;

ft_checkdata(source, 'datatype', 'source')

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source)

cfg.keepindividual = 'yes';
grandavg = ft_sourcegrandaverage(cfg, source, source)

%% format 2

source = [];
source.transform = eye(4);
source.pos = rand(1320,3);
source.pow = rand(1320,1);
source.inside = 1:660;
source.outside = 661:1320;

ft_checkdata(source, 'datatype', 'source')

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source)

cfg.keepindividual = 'yes';
grandavg = ft_sourcegrandaverage(cfg, source, source)

%% format 3

source = [];
source.dim = [10 11 12];
source.transform = eye(4);
source.pow = rand(10*11*12,1);
source.inside = 1:660;
source.outside = 661:1320;

ft_checkdata(source, 'datatype', 'source')

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source)

cfg.keepindividual = 'yes';
grandavg = ft_sourcegrandaverage(cfg, source, source)


