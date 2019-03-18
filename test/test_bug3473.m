function test_bug3473

% WALLTIME 00:10:00
% MEM 2gb

%TEST ft_prepare_mesh 

filename = '/home/common/matlab/fieldtrip/data/ftp/tutorial/epilepsy/case1/ctf_data/case1.pos';

polhemus = ft_read_headshape(filename);
polhemus = ft_convert_units(polhemus,'mm');

template = ft_read_vol('standard_bem.mat');
template = ft_convert_units(template,'mm');

cfg = [];
cfg.template.headshape      = polhemus;
cfg.checksize               = inf;
cfg.individual.headmodel    = template;
cfg                         = ft_interactiverealign(cfg);
template                    = ft_transform_geometry(cfg.m,template);

defaced_template                = template;
cfg                             = [];
defaced_template.bnd(1).unit    = 'mm';
defaced                         =  ft_defacemesh(cfg,defaced_template.bnd(1));

defaced_template.bnd(1).pos = defaced.pos;
defaced_template.bnd(1).tri = defaced.tri;

cfg             = [];
cfg.headshape   = polhemus;
cfg.template    = defaced_template.bnd(1);
cfg.method      = 'fittemplate';
fitted          = ft_prepare_mesh(cfg, template.bnd);

cfg = [];
cfg.method = 'bemcp';
headmodel_bem = ft_prepare_headmodel(cfg, fitted);
