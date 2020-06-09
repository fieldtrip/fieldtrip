function test_bug3473

% WALLTIME 00:30:00
% MEM 2gb
% DEPENDENCY ft_prepare_mesh

addpath(genpath(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3473/cpd')));
filename = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/epilepsy/case1/ctf_data/case1.pos');

polhemus = ft_read_headshape(filename);
polhemus = ft_convert_units(polhemus,'mm');

template = ft_read_vol('standard_bem.mat');
template = ft_convert_units(template,'mm');

m = [0.0000    0.9000         0   40.5000;
    -0.8966    0.0000    0.0784         0;
    0.0784   -0.0000    0.8966   18.0000;
    0         0         0         1.0000];
template = ft_transform_geometry(m, template);

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3473/defaced_polhemus.mat'))
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3473/defaced_template.mat'))

ft_plot_mesh(defaced_template);
ft_plot_mesh(defaced_polhemus);
% Method 1
cfg=[];
cfg.method='singlesphere';
template_sphere = ft_prepare_headmodel(cfg, defaced_template);

cfg=[];
cfg.method = 'singlesphere';
polhemus_sphere = ft_prepare_headmodel(cfg, defaced_polhemus);

scale = polhemus_sphere.r/template_sphere.r;

T2 = [1 0 0 template_sphere.o(1);
    0 1 0 template_sphere.o(2);
    0 0 1 template_sphere.o(3);
    0 0 0 1                ];

T1 = [1 0 0 -template_sphere.o(1);
    0 1 0 -template_sphere.o(2);
    0 0 1 -template_sphere.o(3);
    0 0 0 1                 ];

S  = [scale 0 0 0;
    0 scale 0 0;
    0 0 scale 0;
    0 0 0 1 ];

transformation = T1*S*T2;

template_t_sphere = ft_transform_geometry(transformation, template);

% Method 2
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3473/fiducials.mat'))

cfg             = [];
cfg.method      = 'singlesphere';
template_sphere = ft_prepare_headmodel(cfg, fiducials);

polhemus.fid.pos = polhemus.fid.pos+randn(3,3)/100;

cfg              = [];
cfg.method      = 'singlesphere';
polhemus_sphere = ft_prepare_headmodel(cfg, polhemus.fid);

scale = polhemus_sphere.r/template_sphere.r;

T2 = [1 0 0 template_sphere.o(1);
    0 1 0 template_sphere.o(2);
    0 0 1 template_sphere.o(3);
    0 0 0 1                ];

T1 = [1 0 0 -template_sphere.o(1);
    0 1 0 -template_sphere.o(2);
    0 0 1 -template_sphere.o(3);
    0 0 0 1                 ];

S  = [scale 0 0 0;
    0 scale 0 0;
    0 0 scale 0;
    0 0 0 1 ];

transformation = T1*S*T2;

template_fiducials = ft_transform_geometry(transformation, template);

% Method 3

cfg              = [];
cfg.headshape    = polhemus;
cfg.template     = defaced_template;
cfg.method       = 'fittemplate';
template_surface = ft_prepare_mesh(cfg, template.bnd);

% Volume conduction models

cfg              = [];
cfg.conductivity = [0.33 0.0042 0.33];
cfg.method       = 'openmeeg';
headmodel_sphere = ft_prepare_headmodel(cfg, template_t_sphere.bnd);

cfg               = [];
cfg.conductivity = [0.33 0.0042 0.33];
cfg.method        = 'openmeeg';
headmodel_fiducials = ft_prepare_headmodel(cfg, template_fiducials);

cfg               = [];
cfg.conductivity = [0.33 0.0042 0.33];
cfg.method        = 'openmeeg';
headmodel_surface = ft_prepare_headmodel(cfg, template_surface);

%

cfg                          = [];
cfg.method                   = 'singlesphere';
headmodel_singleshell_sphere = ft_prepare_headmodel(cfg, template_t_sphere.bnd(3));

cfg                          = [];
cfg.method                   = 'singlesphere';
headmodel_singleshell_sphere = ft_prepare_headmodel(cfg, template_fiducials.bnd(3));

cfg                          = [];
cfg.method                   = 'singlesphere';
headmodel_singleshell_sphere = ft_prepare_headmodel(cfg, template_surface(3));
