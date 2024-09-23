function [Gain, cfg] = test_duneuro_bst

ft_hastoolbox('duneuro', 1);
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
segprob.dim = size(segprob.brain);
segprob.unit = 'cm';
segprob.coordsys = 'ctf';
segprob.transform = eye(4);
segprob.transform(1,4) = -0.5;
segprob.transform(2,4) = -0.5;
segprob.transform(3,4) = -0.5;

% hexa mesh
cfg = [];
cfg.shift = 0.3;
cfg.method = 'hexahedral';
mesh_vol_hex = ft_prepare_mesh(cfg, segprob);
mesh_vol_hex = ft_convert_units(mesh_vol_hex,'m');

% tetra mesh
cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);
mesh_vol_tet = ft_convert_units(mesh_vol_tet, 'm');

%% define sensors

% i manually passed 5 coils and fixed projections
% or maybe take some from a ctf file? something like this (from test_pull1377.m)

% for MEG data + sensor info
%
% load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'), 'data');
% datameg = data;
% clear data

coils = [5 5 12; 5 12 5; -2 5 5; 12 5 5; 5 -2 5]./100;
projections = [0 0 1 ; 0 1 0; -1 0 0; 1 0 0; 0 -1 0 ];

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
dip_pos = [5.5 6.5 6.5; 5.5 5.5 3.5; 3.5 5.5 5.5]./100;
dip_mom = [0 1 0; 0 0 -1; -1 0 0 ];


%%%%% Here starts the part that is adjusted from bst_duneuro

cfg = duneuro_defaults;

cfg.modality = sens.type; % FIXME should also accommodate eeg and meeg;

if ~isfield(cfg,'isotrop'); cfg.isotrop     = 1; end                       % important to specify in order to write the correct file format (1 will use MSH, 0 will use Cauchy)
if ~isfield(cfg,'lfAvrgRef'); cfg.lfAvrgRef = 0; end                       %  compute average reference 1, otherwise the electrode 1 is the reference and set to 0

%% ===== PREPARE THE INTERFACE TO DUNEURO =====
%% 0 - Initialisation for duneuro computation (configuration used by the duneuro interface)

usemesh = 'hex';
switch usemesh
  case 'tet'
    mesh = mesh_vol_tet;
    cfg.dnMeshElementType = 'tetrahedron';

  case 'hex'
    mesh = mesh_vol_hex;
% reorder the hex -> FIXME this should be done in the headmodel function

    mesh.hex(:,[3 4 7 8]) = mesh.hex(:,[4 3 8 7]);
    cfg.dnMeshElementType = 'hexahedron';


end

%% ===== DUNEURO =====

cfg.deleteOutputFolder = 0; % brainstorm will manage the rest

cfg.advancedMode = 0;

cond = [0.33, 0.43, 0.53];

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
cfg.sourceSpace       = dip_pos'; %NOTE: transpose is intended, due to how the dipoles have been defined originally in the test_pull1663 function (positions in the columns)

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


cfg.conductivity = cond; % TODO : we can do it from here instead to modify the bst file


cfg.modality = 'meg'; % for now, FIXME


%% 5- Duneuro-Brainstorm interface
%% ===== DUNEURO LEADFIELD COMPUTATION   =====

%% 1- Prepare the head model : 
% Write the head file according to the configuration cfg
filename_headmodel = fullfile(cfg.outputpath, 'headmodel'); % extension is dealt with in the function
filename_headmodel = duneuro_write_headmodel(filename_headmodel, mesh, cfg); % output is needed to get the extension

%% 2- Prepare the Source Model : Same for EEG and MEG 
% Write the source/dipole file
filename_dipoles = fullfile(cfg.outputpath,'dipole.txt');
duneuro_write_sourcemodel(filename_dipoles, cfg.sourceSpace);

%% 3- Prepare the Sensor Model : 
filenames = {fullfile(cfg.outputpath, 'coilpos.txt') fullfile(cfg.outputpath, 'coilori.txt')};
cfg.filename_coilpos = filenames{1};
cfg.filename_coilori = filenames{2};
duneuro_write_sensors(sens, filenames);

%% 4- Prepare the Conductivity/tensor Model
filename_cond = fullfile(cfg.outputpath, 'cond');
filename_cond = duneuro_write_conductivity(filename_cond, 'conductivity', cond);

cfg.filename_headmodel = filename_headmodel;
cfg.filename_dipoles   = filename_dipoles;
cfg.filename_cond      = filename_cond;

%%  5- Prepare the Duneuro Configuration file / the minifile
duneuro_write_minifile(cfg, fullfile(cfg.outputpath, cfg.duneuro_configuration_filename));

if ismac
  appname = 'bst_duneuro_meeg_mac64.app';
elseif isunix
  appname = 'bst_duneuro_meeg_linux64.app';
elseif ispc
  appname = 'bst_duneuro_meeg_win64.exe';
end
appname = which(appname);

%%%% @@ Run Duneuro @@ %%%%%%
system([appname ' ' fullfile(cfg.outputpath, cfg.duneuro_configuration_filename)]);

%% 7- Read the lead field matrix (EEG or/and MEG)
cfg = duneuro_read_leadfield(cfg);

%% 8- Post Process the leadfield
cfg = duneuro_lf(cfg, sens);

%% 6- Read the lead field
% fill the bad channel with nan (not possible within the duneuro)
Gain = cfg.meg.lf;

