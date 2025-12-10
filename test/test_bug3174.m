function test_bug3174

% WALLTIME 00:20:00
% MEM 1gb
% DEPENDENCY ft_volumerealign
% DATA private

load(dccnpath('/project/3031000.02/test/bug3174.mat'));

% co-register the CT to the MRI
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.viewresult  = 'yes'; % view realignment result
cfg.coordsys    = 'tal';
ct2              = ft_volumerealign(cfg, cttal, mri);

