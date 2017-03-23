function test_bug3174

% WALLTIME 00:20:00
% MEM 2gb

% TEST ft_volumerealign

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3174.mat'));

% co-register the CT to the MRI
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.viewresult  = 'yes'; % view realignment result
cfg.coordsys    = 'tal';
ct2              = ft_volumerealign(cfg, cttal, mri);

