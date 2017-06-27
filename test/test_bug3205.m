function test_bug3205

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_electroderealign moveinward

elec = ft_read_sens('standard_1020.elc');

% it would be better to remove '-T3','-T4','-T5','-T6'
% as these overlap with electrodes at the same location with a different name

cfg            = [];
cfg.method     = 'moveinward';
cfg.moveinward = 10;
cfg.elec       = elec;
elec_moved1    = ft_electroderealign(cfg);

cfg            = [];
cfg.method     = 'moveinward';
cfg.moveinward = 10;
elec_moved2    = ft_electroderealign(cfg, elec);

assert(isequaln(rmfield(elec_moved1, 'cfg'), rmfield(elec_moved2, 'cfg')));
