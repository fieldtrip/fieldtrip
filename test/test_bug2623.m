function test_bug2623

% MEM 1500mb
% WALLTIME 00:10:00

% Script tests whether parcellated data can be visualised in ft_sourceplot
% 1. parcellated data is unparcelled
%    parcellated2source is a subfunction inside ft_checkdata
% 2. unparcelled data is visualised using ft_sourceplot

% load data ('source'), and parcellation 
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2623.mat'));

% (1) parcellate
cfg = [];
cfg.method       = 'mean';
cfg.parameter    = 'avg';
dataparcellated  = ft_sourceparcellate(cfg, source, parcellation);

% (2) unparcellate
dataunparcellated = ft_checkdata(dataparcellated, 'datatype','source');

% (3) visualise
cfg = [];
cfg.funparameter = 'avg';
ft_sourceplot(cfg, dataunparcellated); 

