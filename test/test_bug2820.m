% function test_bug2820

% WALLTIME 00:10:00
% MEM 2gb

% it should be is2Dana && ~is2Dfun, so start with a surface
[ftver, ftpath] = ft_version;
surffile = fullfile(ftpath, 'template', 'anatomy', 'surface_pial_both.mat');
anatomical = ft_read_headshape(surffile);

functional     = [];
functional.dim = [256 256 256];
functional.transform = eye(4);
functional.transform(1,4) = -128.5; % place the center of the cube in the origin
functional.transform(2,4) = -128.5;
functional.transform(3,4) = -128.5;
functional.pow = randn(functional.dim);

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg, functional);

% see whether they fit
ft_determine_coordsys(functional, 'interactive', 'no');
ft_plot_mesh(anatomical, 'edgecolor', 'none', 'facecolor', 'skin');
camlight

cfg = [];
cfg.parameter = 'pow';
% cfg.interpmethod  default does not cause  aproblem
interp1 = ft_sourceinterpolate(cfg, functional, anatomical);

cfg = [];
cfg.parameter = 'pow';
cfg.interpmethod = 'project'; % needed to get to the bug
interp2 = ft_sourceinterpolate(cfg, functional, anatomical);

%% the figures should be specled green (on MATLAB 2014b and later)
figure; ft_plot_mesh(interp1, 'vertexcolor', interp1.pow, 'edgecolor', 'none'); camlight; colorbar
figure; ft_plot_mesh(interp2, 'vertexcolor', interp2.pow, 'edgecolor', 'none'); camlight; colorbar

%% use the high-level functionality to get the same figure, this will call ft_sourceinterpolate
cfg = [];
cfg.surffile = surffile;
cfg.funparameter = 'pow';
cfg.method = 'surface';
ft_sourceplot(cfg, functional);

%% another approach to use the high-level functionality to get the same figure, again calling ft_sourceinterpolate internally
cfg = [];
cfg.funparameter = 'pow';
cfg.method = 'surface';
ft_sourceplot(cfg, functional, anatomical);


