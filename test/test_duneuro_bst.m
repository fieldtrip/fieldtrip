function Gain = test_duneuro_bst

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
mesh_vol_hex = ft_convert_units(mesh_vol_hex,'m');

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
hold on
ft_plot_mesh(mesh_vol_hex, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
view(120, 30)


% tetra mesh
cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);
mesh_vol_tet = ft_convert_units(mesh_vol_tet, 'm');

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

coils = [5 5 12; 5 12 5; -2 5 5; 12 5 5; 5 -2 5]./100;
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
dip_pos = [5.5 6.5 6.5; 5.5 5.5 3.5; 3.5 5.5 5.5]./100;
dip_mom = [0 1 0; 0 0 -1; -1 0 0 ];

figure
ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', [5 5 5], 'style', 'intersect');
hold on
quiver3(dip_pos(:,1),dip_pos(:,2),dip_pos(:,3),dip_mom(:,1),dip_mom(:,2),dip_mom(:,3),'bo')
view(-200, 15)

%%%%% Here starts the part that is adjusted from bst_duneuro

% Intialize variables
Dims = 3;
nv = size(dip_pos,1);
Gain = NaN * zeros(size(sens.coilpos,1), Dims * nv);errMsg = '';

% Save current folder
str = which('bst_unique_readme.txt','-all'); % <== later use : bst_fullfile(bst_get('BrainstormUserDir'), 'bst_duneuro');
filepath = fileparts(str{1});                      % < == //                       //                   //                    //

curdir = pwd;
solverMethod = 'duneuro';
% isEeg  = strcmpi(OPTIONS.EEGMethod, solverMethod )  && ~isempty(OPTIONS.iEeg);
% isMeg  = strcmpi(OPTIONS.MEGMethod, solverMethod)  && ~isempty(OPTIONS.iMeg);
% isEcog = strcmpi(OPTIONS.ECOGMethod, solverMethod) && ~isempty(OPTIONS.iEcog);
% isSeeg = strcmpi(OPTIONS.SEEGMethod, solverMethod) && ~isempty(OPTIONS.iSeeg);

% Get temp folder
TmpDir = tempdir;

% Open log file
logFile = fullfile(TmpDir, 'duneuro_log.txt');
fid_log = fopen(logFile, 'w');

%% ===== PREPARE THE INTERFACE TO DUNEURO =====
%% 0 - Initialisation for duneuro computation (configuration used by the duneuro interface)
%  find the bst_duneuro_toolbox path

usemesh = 'tet';
switch usemesh
  case 'tet'
    mesh = mesh_vol_tet;
    numberOfEdges = 4;
  case 'hex'
    mesh = mesh_vol_hex;
    numberOfEdges = 8;
end

% Get number of  layers
numberOfLayer =  numel(mesh.tissuelabel);

% Get mesh element type
if numberOfEdges == 4
    dnMeshElementType = 'tetrahedron';
elseif numberOfEdges == 8
    dnMeshElementType = 'hexahedron';
else
    error('Mesh element type is not defined')
end

% % Get the mesh method generation
% load(OPTIONS.FemFiles, 'Comment');
% if contains(Comment,'iso2mesh')
%     meshMethod = 'iso2mesh';
% elseif contains(Comment,'SimNibs')
%     meshMethod = 'SimNibs';
% else
    meshMethod = 'others';
% end

% Get the default conductivity values
% FIXME OPTIONS.Conductivity  = get_standard_conductivity((numberOfLayer));
% Get the names of the tissues
TissueLabels = mesh.tissuelabel;
OPTIONS.FemNames = TissueLabels;

%% ===== DUNEURO =====

cfg = []; % the most important variable that will pass all the variable to prepare the interface to duneuro
cfg.currentPath = curdir;
cfg.pathOfTempOutPut = TmpDir;
cfg.pathOfDuneuroToolbox = filepath; % USE : bst_fullfile(bst_get('BrainstormUserDir'), 'bst_duneuro');
cfg.displayComment  = 0; % update the bst_progress text
cfg.dnMeshElementType = dnMeshElementType;
cfg.deleteOutputFolder = 0; % brainstorm will manage the rest
cfg.meshMethod = meshMethod;

cfg.advancedMode = 0;

cond = [0.33, 0.43, 0.53];
OPTIONS.FemConductivity = cond;

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

% % Get the modality from the OPTIONS structure
% if isEeg; cfg.modality = 'eeg'; goodChannel = OPTIONS.iEeg; end
% if isMeg; cfg.modality = 'meg'; goodChannel = OPTIONS.iMeg; end
% if isEeg && isMeg; cfg.modality = 'meeg'; goodChannel = [OPTIONS.iMeg OPTIONS.iEeg]; end
% if isEcog; cfg.modality = 'ecog'; goodChannel = OPTIONS.iEcog; end
% if isSeeg; cfg.modality = 'seeg'; goodChannel = OPTIONS.iSeeg; end

%% ===== DUNEURO PREPARE GEOMETRY =====
%% 1 - Head Model

cfg.node = mesh.pos;
try
  cfg.elem = [mesh.tet mesh.tissue];
catch
  cfg.elem = [mesh.hex mesh.tissue];
end
cfg.tissuLabel = mesh.tissuelabel;



% Get here the information about the tissu to keep (similar to the openmeeg panel that check the scalp, inner and outer skull )
% vector of one and zeros, the one mean keep the layer acording to the rank
% of the index on the vector.
% [1 1 0 0 ] ==> there is four layer and here the model will keep the label
% 1 and 2 ===> ALWAYS the  label are from inner to outer  ....
% if isEeg; cfg.eegLayerToKeep = eegLayerToKeep; end % < == situation of 1 1 1 1 1 1 <== keep all, not needed if all 1;
% if isMeg; cfg.megLayerToKeep = megLayerToKeep; end
% if isMeg; cfg.megLayerToKeep = megLayerToKeep; end
% if isEcog; cfg.eCogLayerToKeep =  eCogLayerToKeep; end
% if isSeeg; cfg.seegLayerToKeep =  seegLayerToKeep; end
% temporary option for testing the rest of the code

layerToKeep = [1 1 1];
cfg.layerToKeep = layerToKeep; % or may be this one should be sufficient
cfg.shrinkSourceSpace = 0;

cfg.sourceSpace = dip_pos'; %NOTE: transpose is intended, due to how the dipoles have been defined originally in the test_pull1663 function (positions in the columns)

%% 3- Electrode Position / Coils positon / Orientation
% === EEG ===
isEeg = 0;
if isEeg
    %     EegChannel = [];
    %     for iChan = 1: length(OPTIONS.iEeg)
    %         sChan = OPTIONS.Channel(OPTIONS.iEeg(iChan));
    %         if ~isempty(OPTIONS.Channel(iChan).Loc) %
    %             EegChannel(iChan,:) = OPTIONS.Channel(iChan).Loc;
    %         else
    %             EegChannel(iChan,:) = [];
    %         end
    %     end
    % Channel location :
    eegloc = cat(2, OPTIONS.Channel(OPTIONS.iEeg).Loc)';
    cfg.channelLoc = eegloc;
end
% === MEG ===
% get here the infos about the sensors
isMeg = 1;
if isMeg
  % this looks as if it's integration points, i.e. per coil, and
  % then potentially multiple points, MNE style. In FT-analogy this
  % is a single point estimate, and the tra deals with the
  % integration
  MegChannel = [(1:size(sens.coilpos,1))' sens.coilpos sens.coilori ones(size(sens.coilpos,1),1)];
    % 
    % MegChannel = [];
    % for iChan = 1 : length(OPTIONS.iMeg)
    %     sChan = OPTIONS.Channel(OPTIONS.iMeg(iChan));
    %     for iInteg = 1:size(sChan.Loc, 2)
    % 
    %         MegChannel =  [MegChannel; iChan, sChan.Loc(:,iInteg)', sChan.Orient(:,iInteg)', sChan.Weight(iInteg)];
    %     end
    % end
   cfg.MegChannel = MegChannel;
end

%% 4- Conductivity/tensor
% TODO : Adapt this for anisotrpy // maybe could be done from outside and
% put the tensor on OPTIONS.Conductivity or create new field OPTIONS.ConductivityTensor
cfg.conductivity = OPTIONS.FemConductivity; % TODO : we can do it from here instead to modify the bst file


cfg.modality = 'meg'; % for now, FIXME
 %% 5- Duneuro-Brainstorm interface
%% ===== DUNEURO LEADFIELD COMPUTATION   =====
cfg = bst_duneuro_interface_jm(cfg);
% displayOutput = 1; % 1 yes, 0 no
% if displayOutput == 1; tic; cfg = bst_duneuro_interface(cfg); t_direct = toc; end
% if displayOutput == 0; tic; evalc('cfg = bst_duneuro_interface(cfg);'); t_evalc = toc; end % could be longer ... not really ==> todo check for high resolution

%% 6- Read the lead field
% fil the bad channel with nan (not possible within the duneuro)
Gain = cfg.fem_meg_lf;

%if isMeg;     Gain(OPTIONS.iMeg, :) = cfg.fem_meg_lf; end % the final value ==>  cfg.fem_meg_lf = B = Bp + aBs ? check if the case for bst
%if isEeg;     Gain(OPTIONS.iEeg, :) = cfg.fem_eeg_lf; end

%% ===== CLEANUP =====
% TODO : adapt it to bst
% TODO : Delete intermediary files
% /!\ Warning : the file "transferOut.dat" : should be stored somewhere
% it contains the transfer matrix, this matrix could be re-used if the users change,
% ONLY and ONLY, the source space and/or the source model.
%  If the head geometry (mesh), the electrode position and or the conductivity value are changed,
% the forwar model need to be recomputed again (de not reuse the "transferOut.dat")

% Close log file
if ~isempty(fid_log) && (fid_log >= 0) && ~isempty(fopen(fid_log))
    fclose(fid_log);
end
% Go back to initial folder
cd(curdir);

function cfg = bst_duneuro_interface_jm(cfg)

% BST_DUNEURO_INTERFACE : Writes the arguments from bst to duneuro and run the FEM
%
% USAGE:      cfg = bst_duneuro_interface(cfg)
%
% INPUT:
%     - cfg: structure with the fields:
%              Run the example\demo  in order to have the full liste 
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2019 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Takfarinas MEDANI, December 2019; 


%% 0 - Initialisation 
% Set the default parameters used in this function
if ~isfield(cfg,'modality'); cfg.modality  = 'eeg'; end
if ~isfield(cfg,'runFromBst'); cfg.runFromBst = 1; end                    % Works only if called from brainstorm 
if ~isfield(cfg,'currentPath'); cfg.currentPath = pwd; end                 % This function will cd to a temporary file and then return here (pwd) 
if ~isfield(cfg,'useTransferMatrix'); cfg.useTransferMatrix = 1; end % use the transfer matrix is recommended (choice 0 only for duneuroVersion 1)
if ~isfield(cfg,'BstDuneuroVersion'); cfg.BstDuneuroVersion = 2; end                 % 1 previous with separate files, 2the new version combined eeg and meg and binary + txt output,   
if ~isfield(cfg,'isotrop'); cfg.isotrop =1; end                                       % important to specify in order to write the correct file format (1 will use MSH, 0 will use Cauchy)
if ~isfield(cfg,'lfAvrgRef'); cfg.lfAvrgRef = 0; end                              %  compute average reference 1, otherwise the electrode 1 is the reference and set to 0
                                                                                                          %   in that case the for duneuro the electrod 1 is the reference and set to 0                              
if ~isfield(cfg,'pathOfTempOutPut'); cfg.pathOfTempOutPut =  cfg.currentPath; end  

if cfg.runFromBst == 1
    cfg.lfAvrgRef = 0;
    cfg.brainstormModality = cfg.modality; %
    if ~isfield(cfg,'brainstormOutputFolder') % ==> the output folder could be different from the temporary folder in the case where we want to save the transfer for example
        cfg.brainstormOutputFolder = cfg.pathOfTempOutPut;
    end % should be done from the bst    
     % check the for the user writing right 
    [status, values] = fileattrib(cfg.brainstormOutputFolder);
    if ~values.UserWrite 
        fclose('all');
        if cfg.runFromBst == 1; bst_progress('stop'); end
        error('Duneuro can not write the output data  to this folder : %s. \nPlease change to a path with writing permission', cfg.brainstormOutputFolder);
    end
    
    %%%%% UPDATES these values from here :
    % subpart  [brainstorm]
    % if ~isfield(cfg,'brainstormEegSaveTransfer'); cfg.brainstormEegSaveTransfer = 'false'; end %
    % if ~isfield(cfg,'brainstormMegSaveTransfer'); cfg.brainstormMegSaveTransfer = 'false'; end %
end
                                                                            
                                                                   
cfg.displayComment = 0;
cfg.BstDuneuroVersion = 3;


%% ------------- DUNEURO INTERFACE ------------- %%

% cd(fullfile(cfg.pathOfTempOutPut));
%% 1- Prepare the head model : 
% Write the head file according to the configuration cfg
if ~isfield(cfg,'filename');    cfg.filename = fullfile(cfg.pathOfTempOutPut,'head_model'); end  % Use the default name.
cfg = bst_prepare_head_model(cfg);

%% 2- Prepare the Source Model : Same for EEG and MEG 
% Write the source/dipole file
cfg.dipole_filename  = fullfile(cfg.pathOfTempOutPut,'dipole_model.txt');
write_duneuro_dipole_file(cfg.sourceSpace, cfg.dipole_filename);

%% 3- Prepare the Sensor Model : 
cfg.runFromBst = true; % needed for coil file to be written
cfg = bst_prepare_sensor_model(cfg);

%% 4- Prepare the Conductivity/tensor Model
cfg = bst_prepare_conductivity_model(cfg);

%%  5- Prepare the Duneuro Configuration file / the minifile
cfg = bst_prepare_minifile(cfg);

%% 6- Prepapre & Run the Duneuro Application
% cd(fullfile(cfg.pathOfDuneuroToolbox,'bin'));

% define the command line FIXME this does not work ,because it returns by
% default the .exe win64 executable
cfg = bst_set_duneuro_cmd(cfg);

if ~ispc
  cfg.cmd = strrep(cfg.cmd, 'win64', 'linux64.app');
end

%%%% @@ Run Duneuro @@ %%%%%%
%tic; cfg = bst_run_duneuro_cmd(cfg); cfg.time_fem = toc;
system([cfg.cmd ' ' cfg.mini_filename]);

%% 7- Read the lead field matrix (EEG or/and MEG)
cfg = bst_read_duneuro_leadfield(cfg);

%% 8- Post Process the leadfield
cfg.runFromBst = false; % this overrules the implicit 'tra' multiplication performed
cfg = bst_postprocess_lf(cfg);

%% 9- Remove the temporary folder
if cfg.displayComment ==1;disp(['duneruo >>9 - Clean the folder  : ' (fullfile(cfg.pathOfTempOutPut))]);end
if ~isfield(cfg,'deleteOutputFolder'); cfg.deleteOutputFolder = 0; end
if cfg.deleteOutputFolder == 1
    % TODO or leave bst to do it
    disp(['remove the '  (fullfile(cfg.pathOfTempOutPut)) ' from the hard disc']);
end

%% 10- Go back to the work space
if cfg.displayComment ==1;disp(['duneruo >>10 - Going back to the current folder ' cfg.currentPath ]);end
cd(cfg.currentPath)
if ~isfield(cfg,'writeLogFile'); cfg.writeLogFile = 0; end
if cfg.writeLogFile == 1; diary off; end

