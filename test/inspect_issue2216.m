function inspect_issue2216

% WALLTIME 00:10:00
% DEPENDENCY ft_meshrealign ft_interactiverealign ft_plot_ortho ft_plot_slice
% NODATA

load(dccnpath('issue2216.mat')); % contains mri and headshape

%%
% one required change was to plot all three intersecting meshes without the first two being deleted

f = figure;
[hx, hy, hz] = ft_plot_ortho(mri, 'intersectmesh', headshape, 'style', 'subplot');

%%

f = figure;
[hx, hy, hz] = ft_plot_ortho(mri, 'intersectmesh', headshape, 'style', 'intersect');


%%
% another was to allow the headshape to be transparent (or not) without having to use a global "alpha"

f = figure;
ft_plot_headshape(headshape, 'facealpha', 0.5);

%%
% another was to allow the mri to be transparent (or not) without having to use a global "alpha"

f = figure;
ft_plot_ortho(mri, 'style', 'intersect', 'facealpha', 0.5);

%%
% some changes were in ft_interactiverealign
% the following should NOT show the global alpha

cfg = [];
cfg.template.mri = mri;
cfg.template.mristyle = {'facealpha', 1.0};
cfg.individual.headshape = headshape;
cfg.individual.headshapestyle = {'facealpha', 1, 'material', 'dull'}; 
ft_interactiverealign(cfg)

%%
% some changes in ft_meshrealign, which calls ft_interactiverealign

cfg = [];
cfg.method = 'interactive';
cfg.mri = mri;
ft_meshrealign(cfg, headshape)

%%
% this should now also work

cfg = [];
cfg.method = 'interactive';
cfg.headmodel = [];
cfg.headmodel.r = 92;
cfg.headmodel.o = [0 0 40];
ft_meshrealign(cfg, headshape)


%%
% oh, and this also works

cfg = [];
cfg.intersectmesh = headshape;
ft_sourceplot(cfg, mri)
