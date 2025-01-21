function test_bug2316

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_databrowser
% DATA private

% this function currently serves as a placeholder to reproduce bug2316,
% i.e. it crashes MATLAB (2011a, but possibly other versions as well) on
% a PC running Windows 7. on MacOS it seems fine. 

filename = dccnpath('/project/3031000.02/test/bug2316.mat');
load(filename);

cfg = [];
cfg.layout = 'biosemi128.lay';
cfg.viewmode = 'component';
ft_databrowser(cfg, datacomp);
