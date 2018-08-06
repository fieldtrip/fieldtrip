function test_bug3042

% WALLTIME 00:20:00
% MEM 2gb

% TEST ft_read_headshape ft_read_atlas

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug3042/tess_cortex_pial_low.mat');

%%

original = load(filename);

%%

sourcemodel = ft_read_headshape(filename);
figure
ft_plot_mesh(sourcemodel, 'vertexcolor', sourcemodel.curv)

%%

atlas = ft_read_atlas(filename);
figure
ft_plot_mesh(atlas, 'vertexcolor', atlas.desikan_killiany);
colormap colorcube

%%

cfg = [];
cfg.parcellation = 'destrieux';
parcels = ft_sourceparcellate(cfg, sourcemodel, atlas);


