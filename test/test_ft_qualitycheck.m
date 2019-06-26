function test_ft_qualitycheck

% MEM 1500mb
% WALLTIME 00:30:00

% DEPENDENCY ft_qualitycheck

cfg           = [];
cfg.dataset   = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds'); 
cfg.savemat   = 'no';
cfg.visualize = 'no';
cfg.saveplot  = 'no';
ft_qualitycheck(cfg);
