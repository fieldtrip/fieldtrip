function test_pull1663

% MEM 14gb
% WALLTIME 00:60:00
% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates a set of input-structures to be used for testing
% the DUNEuro MEG forward solution.
% The structure of this script is more or less
% 1. create input data (dull segmentation, hex and tet meshes, dipoles, sensors)
% 2. compute and compare MEG leadfield for hex and tet meshes
%   a. create volume conductor (ft_prepare_headmodel)
%   b. create source grid (ft_prepare_sourcemodel)
%   c. compute leadfield (ft_prepare_leadfield)
%   d. compare hex leadfield with tet leadfield

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. create input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create mesh

% from test_pull1427.m
% dull segmentation
segprob = [];
segprob.brain = false(10,10,10); segprob.brain(4:7,4:7,4:7) = true;
segprob.skull = false(10,10,10); segprob.skull(3:8,3:8,3:8) = true;
segprob.scalp = false(10,10,10); segprob.scalp(2:9,2:9,2:9) = true;
%segprob.air = true(10,10,10);
%segprob.air(2:9,2:9,2:9) = false;
%segprob.brain = false(11,11,11); segprob.brain(4:8,4:8,4:8) = true;
%segprob.skull = false(11,11,11); segprob.skull(3:9,3:9,3:9) = true;
%segprob.scalp = false(11,11,11); segprob.scalp(2:10,2:10,2:10) = true;
segprob.dim = size(segprob.brain);
segprob.unit = 'cm';
segprob.coordsys = 'ctf';
segprob.transform = eye(4);
segprob.transform(1,4) = -0.5;
segprob.transform(2,4) = -0.5;
segprob.transform(3,4) = -0.5;

% visualize the segmentation
% it is more difficult to visualize a probabilistic segmentation than an indexed one
segindx = ft_datatype_segmentation(segprob, 'segmentationstyle', 'indexed');

cfg = [];
cfg.funparameter = 'tissue';
cfg.method = 'ortho';
cfg.location = [5 5 5]; % this is the center of the volume, in this plot it will be rounded off to the nearest voxel
ft_sourceplot(cfg, segindx);

% hexa mesh
cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
mesh_vol_hex = ft_prepare_mesh(cfg, segprob);

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
hold on
ft_plot_mesh(mesh_vol_hex, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)


% tetra mesh
cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);


figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
hold on
ft_plot_mesh(mesh_vol_tet, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)

%% define sensors

% i manually passed 5 coils and fixed projections
% or maybe take some from a ctf file? something like this (from test_pull1377.m)

% for MEG data + sensor info
%
% load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'), 'data');
% datameg = data;
% clear data

coils = [5 5 12; 5 12 5; -2 5 5; 12 5 5; 5 -2 5];
projections = [0 0 1 ; 0 1 0; -1 0 0; 1 0 0; 0 -1 0 ];

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
hold on
quiver3(coils(:,1),coils(:,2),coils(:,3),projections(:,1),projections(:,2),projections(:,3),'bo')

meg_labels = cellstr(strings(1,size(coils,1)));
for i=1:size(coils,1)
  meg_labels(i) = cellstr(strcat('meg',num2str(i)));
end

sens = [];
sens.coilpos = coils;
sens.coilori = projections;
sens.chanpos = coils;
sens.chanori = projections;
sens.label = meg_labels;
sens.type = 'meg';
sens.unit = 'mm';
sens.tra = eye(5);
sens = ft_convert_units(sens,'mm');

%% define dipoles
dip_pos = [5.5 6.5 6.5; 5.5 5.5 3.5; 3.5 5.5 5.5];
dip_mom = [0 1 0; 0 0 -1; -1 0 0 ];

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
hold on
quiver3(dip_pos(:,1),dip_pos(:,2),dip_pos(:,3),dip_mom(:,1),dip_mom(:,2),dip_mom(:,3),'bo')
view(-200, 15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. compute the leadfield
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare headmodel

cfg              = [];
cfg.method       = 'duneuro';
cfg.conductivity = [0.33, 0.43, 0.53];  % vector, conductivity values for tissues: check the order here
vol_duneuro_hex  = ft_prepare_headmodel(cfg, mesh_vol_hex);
vol_duneuro_tet  = ft_prepare_headmodel(cfg, mesh_vol_tet);


%% prepare sourcemodel

% hex
cfg = [];
cfg.sourcemodel.pos = dip_pos';
cfg.sourcemodel.inside = ones(size(dip_pos,2),1);
cfg.grad = sens;
cfg.headmodel = vol_duneuro_hex;
sm_duneuro_hex = ft_prepare_sourcemodel(cfg);

% tet
cfg = [];
cfg.sourcemodel.pos = dip_pos';
cfg.sourcemodel.inside = ones(size(dip_pos,2),1);
cfg.grad = sens;
cfg.headmodel = vol_duneuro_tet;
sm_duneuro_tet = ft_prepare_sourcemodel(cfg);


%% prepare leadfield

% hex
cfg                 = [];
cfg.sourcemodel     = sm_duneuro_hex;
cfg.sourcemodel.mom = dip_mom;
cfg.headmodel       = vol_duneuro_hex;
cfg.grad            = sens;
cfg.reducerank      = 3;
out_hex = ft_prepare_leadfield(cfg);
lf_hex = cell2mat(out_hex.leadfield);

% tet
cfg                 = [];
cfg.sourcemodel     = sm_duneuro_tet;
cfg.sourcemodel.mom = dip_mom;
cfg.headmodel       = vol_duneuro_tet;
cfg.grad            = sens;
cfg.reducerank      = 3;
out_tet = ft_prepare_leadfield(cfg);
lf_tet = cell2mat(out_tet.leadfield);

%% compare leadfields

figure, plot(lf_hex,'r'), hold on, plot(lf_tet,'b')
% set a limit for an error
rel_err_perc = 100*(abs(lf_hex-lf_tet)./norm(lf_hex));
assert(max(max(rel_err_perc))<10)

end % main function

