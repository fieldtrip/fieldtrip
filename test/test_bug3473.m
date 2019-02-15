function test_bug3453
filename = '/home/common/matlab/fieldtrip/data/ftp/tutorial/epilepsy/case1/ctf_data/case1.pos'

polhemus = ft_read_headshape(filename);
% template must be a volume conduction model
template = ft_read_vol('standard_bem.mat');

%%
cfg = [];
% for MEG, polhemus is in MEG coordinates, template is not relevant
cfg.template.headshape      = polhemus;
cfg.individual.headmodel    = template;
cfg.checksize               = inf;
cfg                         = ft_interactiverealign(cfg);
template                    = ft_transform_vol(cfg.m,template);
%%

cfg             = [];
cfg.headshape   = polhemus;
cfg.method      = 'fittemplate';
fitted          = ft_prepare_mesh(cfg, template);


