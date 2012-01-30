function test_1297

% TEST test_1297
% TEST ft_volumesegment

mri_nom = ft_read_mri('/home/common/matlab/fieldtrip/data/test/bug1297/orig-nomask.mgz');

cfg           = [];
cfg.coordsys  = 'spm'; 
cfg.output    = {'skullstrip' 'brain'};
seg2           = ft_volumesegment(cfg, mri_nom);


