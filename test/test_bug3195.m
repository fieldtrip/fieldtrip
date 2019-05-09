function test_bug3195

% MEM 8gb
% WALLTIME 00:20:00
% DEPENDENCY ft_prepare_sourcemodel ft_inside_headmodel

%%
% load the data

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3195'));

load('template_FEM.mat')
load('sourcemodel.mat')

% the system matrix has not been computed yet, hence ft_headmodeltype fails
template_FEM.type = 'simbio';

%%
% take a subset of the cortical sheet vertices

sel = randperm(8004);
sel = sel(1:1000);

cfg = [];
cfg.headmodel = template_FEM;
cfg.sourcemodel.pos = sourcemodel.pos(sel,:);
sourcemodel1 = ft_prepare_sourcemodel(cfg);

figure
ft_plot_mesh(sourcemodel1.pos(sourcemodel1.inside,:))

%%
% construct a coarse 3D grid

cfg = [];
cfg.headmodel = template_FEM;
cfg.resolution = 20;
sourcemodel2 = ft_prepare_sourcemodel(cfg);

figure
ft_plot_mesh(sourcemodel2.pos(sourcemodel2.inside,:))
