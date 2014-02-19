zipFilename = 'H:/common/matlab/fieldtrip/data/ftp/tutorial/Subject01.zip';
system(['unzip ' zipFilename ' ' 'Subject01.mri' ' -d ' pwd]);
mri = ft_read_mri('Subject01.mri');

if isunix
  system(['rm ' 'Subject01.mri']);
elseif ispc
  system(['del ' 'Subject01.mri']);
end

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

