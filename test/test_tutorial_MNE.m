function test_tutorial_MNE

% MEM 4gb
% WALLTIME 00:30:00

% TEST ft_volumereslice ft_volumerealign ft_volumesegment

cd(dccnpath('/home/common/matlab/fieldtrip/data'));
mri = ft_read_mri('Subject01.mri');

cfg        = [];
% cfg.method = 'interactive';
cfg.method = 'fiducial'; % the following voxel coords were determined interactive
cfg.fiducial.nas = [89 60 116];
cfg.fiducial.lpa = [29 145 154];
cfg.fiducial.rpa = [143 142 158];
mri        = ft_volumerealign(cfg, mri);

cfg            = [];
cfg.resolution = 1;
cfg.dim        = [256 256 256];
mri_rs         = ft_volumereslice(cfg, mri);

cfg          = [];
% cfg.method   = 'interactive';
cfg.method = 'fiducial'; % the following voxel coords were determined interactive
cfg.coordsys = 'spm';
cfg.fiducial.ac = [ 153 128 125 ];
cfg.fiducial.pc = [ 101 128 138 ];
cfg.fiducial.xzpoint = [ 138 128 180 ];
mri_mni      = ft_volumerealign(cfg, mri_rs);

cfg           = [];
cfg.output    = {'skullstrip' 'brain'};
mri_seg       = ft_volumesegment(cfg, mri_mni);

