function test_bug3473

% WALLTIME 00:10:00
% MEM 2gb

% DEPENDENCY ft_prepare_mesh 

filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/epilepsy/case1/ctf_data/case1.pos');

polhemus = ft_read_headshape(filename);
polhemus = ft_convert_units(polhemus,'mm');

template = ft_read_vol('standard_bem.mat');
template = ft_convert_units(template,'mm');

cfg             = [];
cfg.headshape   = polhemus;
cfg.template    = template.bnd(1);
cfg.method      = 'fittemplate';
fitted          = ft_prepare_mesh(cfg, template.bnd);

cfg = [];
cfg.method = 'bemcp';
headmodel_bem = ft_prepare_headmodel(cfg, fitted);
