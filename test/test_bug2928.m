function test_bug2928

% WALLTIME 00:20:00
% MEM 4gb

% TEST ft_volumesegment ft_volumerealign

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2928/e889.dcm');

mri = ft_read_mri(filename);

cfg = [];
cfg.coordsys = 'neuromag';
cfg.method = 'fiducial';
cfg.fiducial.nas = [39 135 88];
cfg.fiducial.lpa = [128 168 161];
cfg.fiducial.rpa = [128 164 15];
mri_realigned = ft_volumerealign(cfg, mri);

cfg = [];
mri_segmented = ft_volumesegment(cfg, mri_realigned);  % THIS IS WHERE IT SUPPOSEDLY CRASHES

cfg = [];
cfg.funparameter = 'white';
ft_sourceplot(cfg, mri_segmented);
