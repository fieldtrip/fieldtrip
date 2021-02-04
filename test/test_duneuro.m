%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compares duneuro eeg forward calculation with the one implemented in fieldtrip
% by Sophie Schrader, 16.11.20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('/home/schrader/software/fieldtrip')
ft_defaults

addpath('/home/schrader/software/fieldtrip/forward')
addpath('/home/schrader/software/fieldtrip/external/duneuro')

addpath(genpath('path/to/mex'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this should come from fieldtrip repo or home on donders to run the test

path_data = '.../duneuro-tests/src/';
filename_grid = strcat(path_data, 'test_sphere_tet.msh');
filename_nodes = strcat(path_data, 'test_sphere_tet_nodes.txt');
filename_elements = strcat(path_data, 'test_sphere_tet_elements.txt');
filename_labels = strcat(path_data, 'test_sphere_tet_labels.txt');
filename_tensors = strcat(path_data, 'test_sphere_tet.cond');
filename_coils = strcat(path_data, 'test_sphere_coils.txt');
filename_projections = strcat(path_data, 'test_sphere_projections.txt');
filename_dipoles = strcat(path_data, 'test_sphere_dipoles.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% leadfield without fieldtrip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

forward = 'venant';
% 
% %% create driver object
% cfg = [];
% cfg.type = 'fitted';
% cfg.element_type = 'tetrahedron';
% cfg.solver_type = 'cg';
% cfg.volume_conductor.grid.filename = filename_grid;
% cfg.volume_conductor.tensors.filename = filename_tensors;
% cfg.meg.intorderadd= '0' ;
% cfg.meg.type = 'physical';
% driver = duneuro_meeg(cfg);
% 
% %% read coils
% coils = table2array(readtable(filename_coils))';
% projections = table2array(readtable(filename_projections))';
% driver.set_coils_and_projections(coils, projections);
% 
% %% compute transfer matrix
% cfg = [];
% cfg.solver.reduction = '1e-10';
% cfg.solver.intorderadd= '0' ;
% transfer_matrix = driver.compute_meg_transfer_matrix(cfg);
% 
% %% read dipoles
% dipoles = table2array(readtable(filename_dipoles))';
% 
% %% compute lead field matrix
% cfg = [];
% cfg.post_process = 'true';
% cfg.subtract_mean = 'true';
% 
% if(strcmp(forward,'partial_integration'))
%   cfg.source_model.type = 'partial_integration';
% elseif(strcmp(forward,'venant'))
%   cfg.source_model.type         = 'venant';
%   cfg.source_model.initialization      = 'closest_vertex';
%   cfg.source_model.intorderadd     = '2';
%   cfg.source_model.intorderadd_lb  = '2';
%   cfg.source_model.numberOfMoments = '3';
%   cfg.source_model.referenceLength = '20';
%   cfg.source_model.relaxationFactor = '1e-6';
%   cfg.source_model.restrict        = 'true';
%   cfg.source_model.weightingExponent = '1';
%   cfg.source_model.mixedMoments = 'false';
% end
% 
% lead_field = driver.apply_meg_transfer(transfer_matrix, dipoles, cfg);
% 
% c = [127 127 127];
% mu = 4*pi*1e-4;
% 
% index = repmat(1:size(coils,2),3,1);
% coils_singleprojection = coils(:,index);
% singleprojections = reshape(projections(:),[3 size(coils_singleprojection,2)]);
% Bprimary = compute_B_primary(coils_singleprojection', dipoles', singleprojections');
% 
% Bfull = (mu/(4*pi)) * (Bprimary - lead_field);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% leadfield with fieldtrip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create mesh

% the meshes will be created in fieldtrip

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

%tetra
cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);


%% create the headmodel

% for simbio is like this:
% cfg = [];
% cfg.method ='simbio';
% cfg.conductivity = [0.33, 0.43, 0.53]; % order follows mesh.tissuelabel
% ft_prepare_headmodel(cfg, mesh_vol_tet)

cfg                 = [];
cfg.method          = 'duneuro';
% cfg.grid_filename   = filename_grid; %why i pass also the file? i have the info in the msh struct already, right?
% cfg.tensors_filename= filename_tensors; %can i pass values directly? 
%                                       yes, instead of these 2 lines
%                                       write: 
cfg.conductivity = [0.33, 0.43, 0.53];  % vector, conductivity values for tissues
%                                       check the order here


% this should be all hidden %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% it looks to be optional: try
% cfg.duneuro_settings.type            = 'fitted';
% cfg.duneuro_settings.element_type    = 'tetrahedron';
% cfg.duneuro_settings.solver_type     = 'cg';
% 
% cfg.duneuro_settings.electrodes = 'closest_subentity_center';
% cfg.duneuro_settings.subentities = '3';
% 
% if(strcmp(forward,'partial_integration'))
%   cfg.duneuro_settings.forward         = 'partial_integration';
% elseif(strcmp(forward,'venant')) % this is selected as the way to go
%   cfg.duneuro_settings.forward         = 'venant';
%   cfg.duneuro_settings.initialization  = 'closest_vertex';
%   cfg.duneuro_settings.intorderadd     = '2';
%   cfg.duneuro_settings.intorderadd_lb  = '2';
%   cfg.duneuro_settings.numberOfMoments = '3';
%   cfg.duneuro_settings.referenceLength = '20';
%   cfg.duneuro_settings.relaxationFactor = '1e-6';
%   cfg.duneuro_settings.restrict        = 'true';
%   cfg.duneuro_settings.weightingExponent = '1';
%   cfg.duneuro_settings.mixedMoments = 'false';
% end
% 
% cfg.duneuro_settings.post_process    = 'true';
% cfg.duneuro_settings.subtract_mean   = 'true';
% cfg.duneuro_settings.reduction       = '1e-10';
% cfg.duneuro_settings.intorderadd_meg ='0';
% this should be all hidden %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_duneuro         = ft_prepare_headmodel(cfg, msh);



%% define sensors

% i will manually pass 5 coils and fixed projections or maybe take some
% from a ctf file? something like this (from test_pull1377.m)

% for MEG data + sensor info

% load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'), 'data');
% datameg = data;
% clear data

% but i want to extract only sensor info and not the data at this point

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


%% prepare sourcemodel

% this section has to be checked again

dipoles = table2array(readtable(filename_dipoles))';
dip_pos = dipoles(1:3,:);
dip_mom = dipoles(4:6,:);

cfg = [];
cfg.sourcemodel.pos = dip_pos';
cfg.sourcemodel.inside = ones(size(dipoles,2),1); %check what is considered inside, error?
cfg.grad = sens;
cfg.headmodel = vol_duneuro;
[grid_duneuro] = ft_prepare_sourcemodel(cfg);


%% prepare leadfield
avg_meg = []; %not sure that is necessary

cfg                 = [];
cfg.sourcemodel     = grid_duneuro;
cfg.sourcemodel.mom = dip_mom; %here the dipole moments are passed, ft computes the lf for all 3 directions and then multiplies with moment
cfg.headmodel       = vol_duneuro;
cfg.grad            = sens;
cfg.reducerank      = 3;
[grid_duneuro2]     = ft_prepare_leadfield(cfg, avg_meg);

lf_ft = cell2mat(grid_duneuro2.leadfield);

