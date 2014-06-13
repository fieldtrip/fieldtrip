function test_tutorial_MNE

% MEM 4gb
% WALLTIME 00:30:00

cd(dccnpath('/home/common/matab/fieldtrip/data');
mri = ft_read_mri('Subject01.mri');

cfg        = [];
cfg.method = 'interactive';
mri        = ft_volumerealign(cfg, mri_other);
mri_other = ft_determine_coordsys(mri, 'interactive', 'yes');


cfg            = [];
cfg.resolution = 1;
cfg.dim        = [256 256 256];
mrirs          = ft_volumereslice(cfg, mri);

cfg        = [];
cfg.method = 'interactive';
mri_mni    = ft_volumerealign(cfg, mrirs);

cfg           = [];
cfg.coordsys  = 'spm';
cfg.output    = {'skullstrip' 'brain'};
seg           = ft_volumesegment(cfg, mri_mni);

