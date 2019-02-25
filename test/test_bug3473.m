function test_bug3453
filename = '/home/common/matlab/fieldtrip/data/ftp/tutorial/epilepsy/case1/ctf_data/case1.pos';

polhemus = ft_read_headshape(filename);
polhemus = ft_convert_units(polhemus,'mm');
% template must be a volume conduction model
template = ft_read_vol('standard_bem.mat');
template = ft_convert_units(template,'mm');
%%
cfg = [];
% for MEG, polhemus is in MEG coordinates, template is not relevant
cfg.template.headshape      = polhemus;
cfg.checksize               = inf;
cfg.individual.headmodel    = template;
cfg                           = ft_interactiverealign(cfg);
template                    = ft_transform_vol(cfg.m,template);
%%
defaced_template            = template;
cfg = [];
defaced_template.bnd(1).unit = 'mm';
defaced      =  ft_defacemesh(cfg,defaced_template.bnd(1));

defaced_template.bnd(1).pos = defaced.pos;
defaced_template.bnd(1).tri = defaced.tri;

cfg = [];
defaced_polhemus = ft_defacemesh(cfg,polhemus);
%%
cfg             = [];
cfg.headshape   = defaced_polhemus;
cfg.template    = defaced_template.bnd(1);
cfg.method      = 'fittemplate';
fitted          = ft_prepare_mesh(cfg, template.bnd);
%%
ft_plot_mesh(fitted.bnd)
ft_plot_mesh(defaced_polhemus)
