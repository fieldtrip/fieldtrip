function test_bug3205

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_electroderealign moveinward

elec = ft_read_sens('standard_1020.elc');

cfg = [];
cfg.channel = {'all','-T3','-T4','-T5','-T6'};
elec = ft_selectdata(cfg,elec);

cfg            = [];
cfg.method     = 'moveinward';
cfg.moveinward = 10;
cfg.elec       = elec;
elec_moved     = ft_electroderealign(cfg);

