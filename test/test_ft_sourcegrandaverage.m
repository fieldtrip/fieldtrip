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

ft_checkdata(source, 'datatype', 'source');

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source);

cfg.keepindividual = 'yes';
grandavg = ft_sourcegrandaverage(cfg, source, source);

%% format 2

source = [];
source.transform = eye(4);
source.pos = rand(1320,3);
source.pow = rand(1320,1);
source.inside = 1:660;
source.outside = 661:1320;

ft_checkdata(source, 'datatype', 'source');

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source);

cfg.keepindividual = 'yes';
grandavg = ft_sourcegrandaverage(cfg, source, source);

%% format 3

source = [];
source.dim = [10 11 12];
source.transform = eye(4);
source.pow = rand(10*11*12,1);
source.inside = 1:660;
source.outside = 661:1320;

ft_checkdata(source, 'datatype', 'source');

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source);

cfg.keepindividual = 'yes';
grandavg = ft_sourcegrandaverage(cfg, source, source);


%% format 4

source = [];
source.pos = rand(1320,3);
source.time = 1:25;
source.avg.pow = rand(10*11*12,25);
source.inside = 1:660;
source.outside = 661:1320;

ft_checkdata(source, 'datatype', 'source');

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source);
assert(isfield(grandavg, 'time'), 'time field is missing');

cfg.keepindividual = 'yes';
grandavg = ft_sourcegrandaverage(cfg, source, source);
assert(isfield(grandavg, 'time'), 'time field is missing');

%% format 5

source = [];
source.pos = rand(1320,3);
source.freq = 1:6;
source.time = 1:5;
source.avg.pow = rand(10*11*12,6,5);
source.inside = 1:660;
source.outside = 661:1320;

ft_checkdata(source, 'datatype', 'source');

cfg = [];
cfg.parameter = 'pow';
cfg.keepindividual = 'no';
grandavg = ft_sourcegrandaverage(cfg, source, source);
assert(isfield(grandavg, 'freq'), 'freq field is missing');
assert(isfield(grandavg, 'time'), 'time field is missing');

cfg.keepindividual = 'yes';
grandavg = ft_sourcegrandaverage(cfg, source, source);
assert(isfield(grandavg, 'freq'), 'freq field is missing');
assert(isfield(grandavg, 'time'), 'time field is missing');



