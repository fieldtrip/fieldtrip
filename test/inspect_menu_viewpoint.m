function inspect_menu_viewpoint

% Some of the FT_PLOT_XXX functions that return a 3D object support a
% right-mouse-click context menu with which you can select
% top/bottom/left/right/front/back. This functionality requires that the object being
% plotted has a known coordinate system.

% note that the objects don't really go together, so it does not fit 100%

figure

%%

elec = ft_read_sens('GSN-HydroCel-128.sfp');
elec = ft_convert_units(elec, 'mm');
elec.coordsys = 'ras';

ft_plot_sens(elec, 'label', 'label'); axis off

%%

mesh = ft_read_headshape('cortex_8196.surf.gii');
mesh.coordsys = 'mni';

ft_plot_mesh(mesh, 'facealpha', 0.5, 'unit', 'mm')
lighting gouraud

%%
% dipole pointing to anterior and superior

ft_plot_dipole([40 0 50], [0 1 1], 'coordsys', 'ras', 'unit', 'mm', 'axes', true)

%%

ft_plot_topo3d(elec.elecpos, elec.elecpos(:,3), 'coordsys', 'ras', 'unit', 'mm', 'axes', true)

