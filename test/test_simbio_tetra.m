function test_simbio_tetra
%
% MEM 6gb
% WALLTIME 1:00:00
% DEPENDENCY ft_prepare_mesh ft_prepare_headmodel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates an hexahedral and tetrahedral volumetric mesh from a 
% two-compartments volume to be passed as input to ft_prepare_headmodel with 
% the method 'simbio'.

segm = [];
segm.brain = logical(floor(3*rand(10,10,10)));
segm.scalp = ~segm.brain;
segm.dim = size(segm.brain);
segm.unit = 'mm';
segm.coordsys = 'ctf';
segm.transform = eye(4);


cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
cfg.resolution = 1; 
mesh_vol_hex = ft_prepare_mesh(cfg, segm);

ft_plot_mesh(mesh_vol_hex, 'surfaceonly', false);
view(120, 30)

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [0.33, 0.43]; % order follows mesh.tissuelabel
ft_prepare_headmodel(cfg, mesh_vol_hex)

%%

cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segm);

figure
ft_plot_mesh(mesh_vol_tet, 'surfaceonly', false);
view(120, 30)

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [0.33, 0.43]; % order follows mesh.tissuelabel
ft_prepare_headmodel(cfg, mesh_vol_tet)

%%%% ERROR: 
% Error using ft_headmodel_simbio (line 65)
% No element indices declared!
end