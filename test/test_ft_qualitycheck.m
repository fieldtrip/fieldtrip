function test_ft_qualitycheck

% TEST test_ft_qualitycheck ft_qualitycheck

pwdir = pwd;
cd('/home/common/matlab/fieldtrip/data/test/original/meg/');

cd('ctf275');
cfg           = [];
cfg.dataset   = 'A0132_Aud-Obj-Recognition_20051115_02.ds'; 
cfg.savemat   = 'no';
cfg.visualize = 'no';
cfg.saveplot  = 'no';
ft_qualitycheck(cfg);

cd(pwdir);
