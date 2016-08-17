function inspect_bug3141

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_defacemesh ft_defacevolume

%% anatomical mri

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));

cfg = [];
defaced = ft_defacevolume(cfg, mri);

cfg = [];
ft_sourceplot(cfg, defaced);


%% head shape

headshape = ft_read_headshape(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.shape'));

cfg = [];
defaced = ft_defacemesh(cfg, headshape);

figure
ft_plot_mesh(defaced);

%% 3D grid source model

% this MATLAB file contains the variable sourcemodel
load(dccnpath('/home/common/matlab/fieldtrip/template/sourcemodel/standard_sourcemodel3d4mm.mat'));

cfg = [];
defaced = ft_defacemesh(cfg, sourcemodel);

figure
ft_plot_mesh(defaced.pos(defaced.inside,:));

%% cortical sheet source model

sourcemodel = ft_read_headshape(dccnpath('/home/common/matlab/fieldtrip/template/sourcemodel/cortex_8196.surf.gii'));

cfg = [];
defaced = ft_defacemesh(cfg, sourcemodel);

figure
ft_plot_mesh(defaced);
camlight
lighting phong
