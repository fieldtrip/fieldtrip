function test_ft_meshrealign

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_meshrealign ft_read_mri ft_read_sens ft_prepare_mesh
% DATA no

[ftver, ftpath] = ft_version;
templatedir  = fullfile(ftpath, 'template');

mrifilename   = fullfile(templatedir, 'headmodel', 'standard_seg.mat');
elecfilename  = fullfile(templatedir, 'electrode', 'standard_1020.elc');

mri = ft_read_mri(mrifilename);
cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
cfg.downsample = 10;
mesh = ft_prepare_mesh(cfg,mri);

elec = ft_read_sens(elecfilename);

cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.fiducial.nas = elec.elecpos(3,:);
cfg.fiducial.lpa = elec.elecpos(1,:);
cfg.fiducial.rpa = elec.elecpos(2,:);
meshout = ft_meshrealign(cfg, mesh);