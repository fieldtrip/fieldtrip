function test_bug3110

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY
% DATA private

% this functions tests the behavior of ft_redefinetrial with comp input

cfg = [];
cfg.inputfile = dccnpath('/project/3031000.02/test/latest/raw/meg/preproc_ctf275.mat');
cfg.method    = 'pca';
comp = ft_componentanalysis(cfg);


cfg = [];
cfg.length = 0.5;
cfg.overlap = 0.5;
comp2 = ft_redefinetrial(cfg, comp);

assert(isfield(comp2, 'topo'));

