function test_bug3042

% WALLTIME 00:20:00
% MEM 2gb

% TEST test_bug3042
% TEST ft_read_headshape ft_read_atlas

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3042/tess_cortex_pial_low.mat');

%%
original = load(filename);

%%
sourcemodel = ft_read_headshape(filename);
figure
ft_plot_mesh(sourcemodel, 'vertexcolor', sourcemodel.curv)

%%
parcellation = ft_read_atlas(filename);
figure
ft_plot_mesh(parcellation, 'vertexcolor', parcellation.desikan_killiany);
colormap colorcube



