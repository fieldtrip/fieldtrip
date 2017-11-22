function failed_bug3177

% WALLTIME 0:20:00
% MEM 3gb

% TEST test_bug3177
% TEST ft_electroderealign mesh2edge poly2tri

%%

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/headmodel_fem/vol.mat'))


%% project electrodes
% this converts on the fly the hexaehders into a polygonal surface mesh and subsequently into a triangulated surface mesh

cfg                 = [];
cfg.elec            = ft_read_sens('standard_1020.elc');
cfg.headshape       = vol;
cfg.method          = 'project';
elec_proj           = ft_electroderealign(cfg);
