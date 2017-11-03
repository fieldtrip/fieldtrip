filename = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri');
mri = ft_read_mri(filename);

ft_file = which('ft_defaults');
[p,f,e] = fileparts(ft_file);

cfg = [];
cfg.spmversion = 'spm2';
n2 = ft_volumenormalise(cfg, mri);

rmpath(fullfile(p,'external','spm2'));
cfg.spmversion = 'spm8';
n8 = ft_volumenormalise(cfg, mri);

rmpath(fullfile(p,'external','spm8'));
cfg.spmversion = 'spm12';
n12 = ft_volumenormalise(cfg, mri);

cfg.spmmethod = 'new';
n12new = ft_volumenormalise(cfg, mri);
