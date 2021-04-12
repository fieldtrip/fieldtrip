try

load sourcemodel.mat % this contains "tmp"

npos = 2472;
functional = keepfields(tmp, {'pos', 'inside', 'unit'});
functional.pow = randn(npos, 1);

anatomical = ft_read_mri('single_subj_T1_1mm.nii');

%%

cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'sphere_avg';
cfg.sphereradius = 2;
interp = ft_sourceinterpolate(cfg, functional, anatomical);

%%

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg, interp);

end
