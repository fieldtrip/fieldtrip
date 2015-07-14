function test_bug2928

% WALLTIME 00:20:00
% MEM 4gb

% TEST test_bug2928
% TEST ft_volumesegment ft_volumerealign

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2928/e889.dcm');

mi = ft_read_mri(filename);

cfg = [];
cfg.coordsys = 'neuromag';
cfg.method = 'fiducial';
cfg.fidicial.nas = [39 135 88];
cfg.fidicial.lpa = [128 168 161];
cfg.fidicial.rpa = [128 164 15];
mri_realigned = ft_volumerealign(cfg, mri);

cfg = [];
seg = ft_volumesegment(cfg, mri_realigned);  % THIS IS WHERE IT SUPPOSEDLY CRASHES

cfg = [];
cfg.funparameter = 'white';
ft_sourceplot(cfg, seg);
