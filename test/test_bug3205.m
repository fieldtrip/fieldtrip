function test_bug3205

% WALLTIME 00:00:05
% MEM 5mb
% TEST test_bug3205
% TEST ft_electroderealign

elec = ft_read_sens('standard_1020.elc');
cfg.channel = {'all','-T3','-T4','-T5','-T6'};
elec = ft_selectdata(cfg,elec);

cfg            = [];
cfg.method     = 'moveinward';
cfg.moveinward = 10;
cfg.elec       = elec;
elec_moved     = ft_electroderealign(cfg);

