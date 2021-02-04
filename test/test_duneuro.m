function test_duneuro

% MEM 12gb %%%check
% WALLTIME 1:30:00 %%%check
% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates a set of input-structures to be used for testing
% the duneuro MEG forward solution.
% The structure of this script is more or less
% 1. create input data (dull segmentation, hex and tet meshes, dipoles, sensors)
% 2. compute MEG leadfield for hex and tet meshes
%   a. create volume conductor (ft_prepare_headmodel) 
%   b. create source grid (ft_prepare_sourcemodel)
%   c. compute leadfield (ft_prepare_leadfield)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In general, variables in one function workspace are not available to other
% functions. However, nested functions can access and modify variables in the
% workspaces of the functions that contain them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adapted to be ft compliant from script by Sophie Schrader, 16.11.20

addpath /home/neurophys/marpia/fieldtrip/
ft_defaults

% prevent errors from cfg.mne.keepleadfield, etc
global ft_default
ft_default.checkconfig = 'loose';

% addpath(('/home/common/matlab/fieldtrip/external/duneuro/'))
addpath(genpath('/home/neurophys/marpia/fieldtrip/external/duneuro')) %should be moved somewhere
% addpath /home/neurophys/marpia/Dune2.6/build-release6/duneuro-matlab/src

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
% segprob.air = true(10,10,10);
% segprob.air(2:9,2:9,2:9) = false;
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

% visualize the segmentation %%%%%%%%%%%%%%%%%%%%

% it is more difficult to visualize a probabilistic segmentation than an indexed one
segindx = ft_datatype_segmentation(segprob, 'segmentationstyle', 'indexed');

cfg = [];
cfg.funparameter = 'seg';
cfg.method = 'ortho';
cfg.location = [5 5 5]; % this is the center of the volume, in this plot it will be rounded off to the nearest voxel
ft_sourceplot(cfg, segindx);

% determine the range of the bounding box
[X, Y, Z] = ndgrid(1:segprob.dim(1), 1:segprob.dim(2), 1:segprob.dim(3));
voxpos = ft_warp_apply(segprob.transform, [X(:) Y(:) Z(:)]);
minmaxpos(:,1) = min(voxpos) - 0.5;
minmaxpos(:,2) = max(voxpos) + 0.5;

% visualize the segmentation %%%%%%%%%%%%%%%%%%%%

% hexa
cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
mesh_vol_hex = ft_prepare_mesh(cfg, segprob);

% tetra
cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);

% plot the meshes?

%% define sensors

% i will manually pass 5 coils and fixed projections or maybe take some
% from a ctf file? something like this (from test_pull1377.m)

% for MEG data + sensor info

% load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'), 'data');
% datameg = data;
% clear data

% but i want to extract only sensor info and not the data at this point

% coils = [8 5 -2; 12 5 -2; 8 5 2; 12 5 2; 5 12 5];

coils = table2array(readtable(filename_coils))';
projections = table2array(readtable(filename_projections))';

% why this? fieldtrip wants only one projection per coil, i guess
index = repmat(1:size(coils,2),3,1);
coils_singleprojection = coils(:,index);
singleprojections = reshape(projections(:),[3 size(coils_singleprojection,2)]);

meg_labels = cellstr(strings(1,size(coils_singleprojection,2)));
for i=1:size(coils_singleprojection,2)
  meg_labels(i) = cellstr(strcat('meg',num2str(i)));
end

sens = [];
sens.coilpos = coils_singleprojection';
sens.coilori = singleprojections';
sens.chanpos = coils_singleprojection';
sens.chanori = singleprojections';
sens.label = meg_labels;
sens.type = 'meg';
sens.unit = 'mm';
sens = ft_convert_units(sens,'mm');

%% define dipoles

dipoles = table2array(readtable(filename_dipoles))';
dip_pos = dipoles(1:3,:);
dip_mom = dipoles(4:6,:);

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

% this section has to be checked again

cfg = [];
cfg.sourcemodel.pos = dip_pos';
cfg.sourcemodel.inside = ones(size(dipoles,2),1); %check what is considered inside, error?
cfg.grad = sens;
cfg.headmodel = vol_duneuro_hex;
[sm_duneuro_hex] = ft_prepare_sourcemodel(cfg);

cfg = [];
cfg.sourcemodel.pos = dip_pos';
cfg.sourcemodel.inside = ones(size(dipoles,2),1); %check what is considered inside, error?
cfg.grad = sens;
cfg.headmodel = vol_duneuro_tet;
[sm_duneuro_tet] = ft_prepare_sourcemodel(cfg);


%% prepare leadfield
avg_meg = []; %not sure that is necessary

cfg                 = [];
cfg.sourcemodel     = sm_duneuro_hex;
cfg.sourcemodel.mom = dip_mom; %here the dipole moments are passed, ft computes the lf for all 3 directions and then multiplies with moment
cfg.headmodel       = vol_duneuro_hex;
cfg.grad            = sens;
cfg.reducerank      = 3;
ft_prepare_leadfield(cfg, avg_meg);

cfg                 = [];
cfg.sourcemodel     = sm_duneuro_tet;
cfg.sourcemodel.mom = dip_mom; %here the dipole moments are passed, ft computes the lf for all 3 directions and then multiplies with moment
cfg.headmodel       = vol_duneuro_tet;
cfg.grad            = sens;
cfg.reducerank      = 3;
ft_prepare_leadfield(cfg, avg_meg);

end % main function

