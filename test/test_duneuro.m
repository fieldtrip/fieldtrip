function [lf_bst, cfg] = test_duneuro_bst(nvoxbox)

% MEM 12gb
% WALLTIME 01:00:00
% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis
% DATA private

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
% 3. create the leadfields based on a singlesphere model for reference
% (ballpark numbers should coincide)

ft_hastoolbox('duneuro', 1);

if nargin<1
  nvoxbox = 11;
else
  % having nvoxbox slightly bigger allows for potentially more
  % realistically sized computations
end

if nvoxbox<10
  ft_error('the length of the volume conductor cube should be at least 10');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. create input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create mesh

brainlim = ceil(0.4.*nvoxbox):floor(0.7.*nvoxbox);
skulllim = ceil(0.3.*nvoxbox):floor(0.8.*nvoxbox);
scalplim = ceil(0.2.*nvoxbox):floor(0.9.*nvoxbox);

% from test_pull1427.m
% dull segmentation
segprob = [];
segprob.brain = false(nvoxbox,nvoxbox,nvoxbox); segprob.brain(brainlim,brainlim,brainlim) = true;
segprob.skull = false(nvoxbox,nvoxbox,nvoxbox); segprob.skull(skulllim,skulllim,skulllim) = true;
segprob.scalp = false(nvoxbox,nvoxbox,nvoxbox); segprob.scalp(scalplim,scalplim,scalplim) = true;
segprob.dim = size(segprob.brain);
segprob.unit = 'cm';
segprob.coordsys = 'ctf';
segprob.transform = eye(4);
segprob.transform(1,4) = -0.5;
segprob.transform(2,4) = -0.5;
segprob.transform(3,4) = -0.5;

center  = [nvoxbox nvoxbox nvoxbox]./2;

% visualize the segmentation
% it is more difficult to visualize a probabilistic segmentation than an indexed one
segindx = ft_datatype_segmentation(segprob, 'segmentationstyle', 'indexed');

cfg = [];
cfg.funparameter = 'tissue';
cfg.method   = 'ortho';
cfg.location = center; % this is the center of the volume, in this plot it will be rounded off to the nearest voxel
ft_sourceplot(cfg, segindx);

% hexa mesh
cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
mesh_vol_hex = ft_prepare_mesh(cfg, segprob);
mesh_vol_hex = ft_convert_units(mesh_vol_hex, 'm');

% surface mesh
cfg = [];
cfg.tissue = 'brain';
mesh_surf = ft_prepare_mesh(cfg, segprob);
mesh_surf = ft_convert_units(mesh_surf, 'm');

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', center, 'style', 'intersect');
hold on
ft_plot_mesh(mesh_vol_hex, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)


% tetra mesh
cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);
mesh_vol_tet = ft_convert_units(mesh_vol_tet, 'm');

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', center, 'style', 'intersect');
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

coils = [center(1) center(2) nvoxbox+2; center(1) nvoxbox+2 center(3); -2 center(2) center(3); nvoxbox+2 center(2) center(3); center(1) -2 center(3)]; % 2 cm out of the box
projections = [0 0 1 ; 0 1 0; -1 0 0; 1 0 0; 0 -1 0 ];

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', center, 'style', 'intersect');
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
sens.unit = 'cm';
sens.tra = eye(5);
sens = ft_convert_units(sens,'m');
sens = ft_datatype_sens(sens);
sens.chanunit = repmat({'T'}, [5 1]); % avoid bookkeeping slowdown
sens.chantype = repmat({'megmag'}, [5 1]);

%% define dipoles
dip_pos = [5.5 6.5 6.5; 5.5 5.5 3.5; 3.5 5.5 5.5]./100; % in m
dip_mom = [0 1 0; 0 0 -1; -1 0 0 ];

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', center, 'style', 'intersect');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. compute the leadfield with singlesphere model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare headmodel

cfg              = [];
cfg.method       = 'singlesphere';
cfg.tissue       = 'brain';
vol_ssph          = ft_prepare_headmodel(cfg, segprob);

%% prepare sourcemodel

% hex
cfg = [];
cfg.sourcemodel.pos = dip_pos';
cfg.sourcemodel.inside = ones(size(dip_pos,2),1);
cfg.grad = sens;
cfg.headmodel = vol_ssph;
sm_ssph = ft_prepare_sourcemodel(cfg);

%% prepare leadfield

% hex
cfg                 = [];
cfg.sourcemodel     = sm_ssph;
cfg.sourcemodel.mom = dip_mom;
cfg.headmodel       = vol_ssph;
cfg.grad            = sens;
cfg.reducerank      = 3;
out_ssph = ft_prepare_leadfield(cfg);
lf_ssph = cell2mat(out_ssph.leadfield);

figure, plot(lf_hex,'r'), hold on, plot(lf_tet,'b'), plot(lf_ssph, 'g');

% set a limit for an error
rel_err_perc = 100*(abs(lf_hex-lf_tet)./norm(lf_hex));
assert(max(max(rel_err_perc))<10)

rel_err_perc2 = 100*(abs(lf_ssph-lf_tet)./norm(lf_ssph));
%assert(max(max(rel_err_perc2))<10)



%%%%% Here starts the part that is adjusted from bst_duneuro

cfg = duneuro_defaults;

cfg.modality = 'meg';%sens.type; % FIXME should also accommodate eeg and meeg;

if ~isfield(cfg,'lfAvrgRef'); cfg.lfAvrgRef = 0; end                       %  compute average reference 1, otherwise the electrode 1 is the reference and set to 0

%% ===== PREPARE THE INTERFACE TO DUNEURO =====
%% 0 - Initialisation for duneuro computation (configuration used by the duneuro interface)

%% ===== DUNEURO =====
% % TODO : The conductivities values in the case of the
% % combined model should be the same for eeg and meg
% %% Ask for the layers to keep for the MEG computation
% % We can not use different volme conductor in the case of combined
% % eeg/meg ... otherwise we need to call separately the binaries ....
% % TODO : discuss this point with th the team.
% if strcmpi(OPTIONS.EEGMethod, 'duneuro') && strcmpi(OPTIONS.MEGMethod, 'duneuro')
%     [res, isCancel] = java_dialog('checkbox', ...
%         '<HTML>Select the layers to consider for the combined MEG/EEG head modeling <BR>', 'Select Volume', [], ...
%         TissueLabels, [ones(1, length(TissueLabels))]);
%     if isCancel
%         return
%     end
%     layerToKeep =  res;
% elseif  strcmpi(OPTIONS.EEGMethod, 'duneuro')
%     [res, isCancel] = java_dialog('checkbox', ...
%         '<HTML>Select the layers to consider for the EEG head modeling <BR>', 'Select Volume', [], ...
%         TissueLabels, [ones(1, length(TissueLabels))]);
%     layerToKeep =  res;
%     if isCancel
%         return
%     end
% elseif strcmpi(OPTIONS.MEGMethod, 'duneuro')
%     [res, isCancel] = java_dialog('checkbox', ...
%         '<HTML>Select the layers to consider for the MEG head modeling <BR>', 'Select Volume', [], ...
%         TissueLabels, [1 zeros(1, length(TissueLabels)-1)]);
%     if isCancel
%         return
%     end
%     layerToKeep =  res;
% else
%     % TODO : add similar option in the case of iEEG and sEEG
% end
% 


%% ===== DUNEURO PREPARE GEOMETRY =====
%% 1 - Head Model


% cfg.layerToKeep       = [1 1 1]; % or may be this one should be sufficient
cfg.shrinkSourceSpace = 0;

%% 3- Electrode Position / Coils positon / Orientation
% === EEG ===
% isEeg = 0;
% if isEeg
%     %     EegChannel = [];
%     %     for iChan = 1: length(OPTIONS.iEeg)
%     %         sChan = OPTIONS.Channel(OPTIONS.iEeg(iChan));
%     %         if ~isempty(OPTIONS.Channel(iChan).Loc) %
%     %             EegChannel(iChan,:) = OPTIONS.Channel(iChan).Loc;
%     %         else
%     %             EegChannel(iChan,:) = [];
%     %         end
%     %     end
%     % Channel location :
%     eegloc = cat(2, OPTIONS.Channel(OPTIONS.iEeg).Loc)';
%     cfg.channelLoc = eegloc;
% end
% % === MEG ===
% get here the infos about the sensors
% isMeg = 1;
% if isMeg
%   % this looks as if it's integration points, i.e. per coil, and
%   % then potentially multiple points, MNE style. In FT-analogy this
%   % is a single point estimate, and the tra deals with the
%   % integration
%   MegChannel = [(1:size(sens.coilpos,1))' sens.coilpos sens.coilori ones(size(sens.coilpos,1),1)];
%     % 
%     % MegChannel = [];
%     % for iChan = 1 : length(OPTIONS.iMeg)
%     %     sChan = OPTIONS.Channel(OPTIONS.iMeg(iChan));
%     %     for iInteg = 1:size(sChan.Loc, 2)
%     % 
%     %         MegChannel =  [MegChannel; iChan, sChan.Loc(:,iInteg)', sChan.Orient(:,iInteg)', sChan.Weight(iInteg)];
%     %     end
%     % end
%    cfg.MegChannel = MegChannel;
% end

%% 4- Conductivity/tensor
% TODO : Adapt this for anisotrpy // maybe could be done from outside and
% put the tensor on OPTIONS.Conductivity or create new field OPTIONS.ConductivityTensor

cfg = [];
cfg.duneuro.application = '/home/mrphys/jansch/matlab/fieldtrip/external/duneuro/bst_duneuro_meeg_linux64.app';
cfg.method       = 'duneuro';
cfg.conductivity = [0.33, 0.43, 0.53];  % vector, conductivity values for tissues: check the order here
vol_duneuro2_hex  = ft_prepare_headmodel(cfg, mesh_vol_hex);

cfg                 = [];
cfg.sourcemodel     = sm_duneuro_hex;
cfg.sourcemodel.mom = dip_mom;
cfg.headmodel       = vol_duneuro2_hex;
cfg.grad            = sens;
cfg.reducerank      = 3;
out_hex2 = ft_prepare_leadfield(cfg);
lf_hex2  = cell2mat(out_hex.leadfield);

% set a limit for an error
rel_err_perc = 100*(abs(lf_hex-lf_hex2)./norm(lf_hex));
assert(max(max(rel_err_perc))<10)


