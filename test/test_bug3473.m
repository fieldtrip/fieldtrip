function test_bug3453
filename = 'case1/ctf_data/case1.pos'

polhemus = ft_read_headshape(filename);

template = ft_read_vol('standard_bem.mat');

%%
cfg = [];
% for MEG, polhemus is in MEG coordinates, template is not relevant
 cfg.template.headshape      = polhemus;
 cfg.individual.headmodel    = template;
 cfg.checksize = inf;
[cfg] = ft_interactiverealign(cfg);

template = ft_transform_vol(cfg.m,template);
%%

cfg = [];
cfg.headshape = polhemus;
cfg.method = 'fittemplate';
fitted = ft_prepare_mesh(cfg, template);
% Units have to be consistend
% Conductivity has to be kept
% Stiffnes matrix has to be removed
% template mush be a mesh
% what about a FEM template model?