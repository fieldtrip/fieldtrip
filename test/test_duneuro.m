function test_duneuro(nvoxbox, doplot)

% MEM 12gb
% WALLTIME 01:00:00
% DEPENDENCY ft_prepare_sourcemodel ft_prepare_leadfield
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

if nargin<1 || isempty(nvoxbox)
  nvoxbox = 11;
else
  % having nvoxbox slightly bigger allows for potentially more
  % realistically sized computations
end

if nargin<2
  doplot = false;
end

if nvoxbox<11
  ft_error('the length of the volume conductor cube should be at least 11');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. create input data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create mesh


% from test_pull1427.m
% dull segmentation
segprob = [];
segprob.brain = false(nvoxbox,nvoxbox,nvoxbox); 
segprob.skull = false(nvoxbox,nvoxbox,nvoxbox); 
segprob.scalp = false(nvoxbox,nvoxbox,nvoxbox); 
segprob.dim = size(segprob.brain);
segprob.unit = 'cm';
segprob.coordsys = 'ctf';
segprob.transform = diag([ (10./(nvoxbox-1)).*ones(1,3) 1] );
segprob.transform(1,4) = -5.5;%-nvoxbox./2-0.5;
segprob.transform(2,4) = -5.5;%-nvoxbox./2-0.5;
segprob.transform(3,4) = -5.5;%-nvoxbox./2-0.5;

center  = [0 0 0];

[x,y,z] = ndgrid(1:nvoxbox, 1:nvoxbox, 1:nvoxbox);
pos = ft_warp_apply(segprob.transform, [x(:) y(:) z(:)]);
R   = sqrt(sum(pos.^2,2));
clear x y z

r = (max(pos(:))-min(pos(:)))./2;

segprob.brain(R<=0.5.*r) = true;
segprob.skull(R<=0.7.*r) = true;
segprob.scalp(R<=0.9.*r) = true;

% visualize the segmentation
% it is more difficult to visualize a probabilistic segmentation than an indexed one
segindx = ft_datatype_segmentation(segprob, 'segmentationstyle', 'indexed');
segindx = ft_convert_units(segindx, 'm');

if doplot
  cfg = [];
  cfg.funparameter = 'tissue';
  cfg.method   = 'ortho';
  cfg.location = center; % this is the center of the volume, in this plot it will be rounded off to the nearest voxel
  ft_sourceplot(cfg, segindx);
end

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

if doplot
  figure
  ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', center, 'style', 'intersect');
  hold on
  ft_plot_mesh(mesh_vol_hex, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
  view(120, 30)
end

% tetra mesh
cfg = [];
cfg.method = 'tetrahedral';
mesh_vol_tet = ft_prepare_mesh(cfg, segprob);
mesh_vol_tet = ft_convert_units(mesh_vol_tet, 'm');

if doplot
  figure
  ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', center, 'style', 'intersect');
  hold on
  ft_plot_mesh(mesh_vol_tet, 'surfaceonly', false, 'facecolor', 'none', 'edgecolor', 'm');
  view(120, 30)
end

%% define sensors

% MEG sensor array
coils = [0 0 r+2; 0 r+2 0; -(r+2) 0 0; r+2 0 0; 0 -(r+2) 0]; % 2 cm out of the box
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
sens.unit = 'cm';
sens.tra = eye(5);
sens = ft_convert_units(sens,'m');
sens = ft_datatype_sens(sens);
sens.chanunit = repmat({'T'}, [5 1]); % avoid bookkeeping slowdown
sens.chantype = repmat({'megmag'}, [5 1]);

if doplot
  figure
  ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', center, 'style', 'intersect');
  hold on
  ft_plot_sens(sens);
end

% EEG sensor array
elec = [];
elec.elecpos = [0 0 r; 0 r 0; -r 0 0; r 0 0; 0 -r 0].*0.9;
elec.type    = 'eeg';
elec.unit    = 'cm';
elec.label   = {'eeg1';'eeg2';'eeg3';'eeg4';'eeg5'};
elec.chanunit = {'V' 'V' 'V' 'V' 'V'}';
elec.chantype = {'eeg' 'eeg' 'eeg' 'eeg' 'eeg'}';
elec = ft_convert_units(elec, 'm');
elec = ft_datatype_sens(elec);


%% define dipoles

% -> just picking randomly yields very poor solutions (without tweaking)
% when using hexagonal mesh
%dip_pos = [rand(1,3); rand(1,3).*[1 1 -1]; rand(1,3).*[-1 1 1]];
%dip_pos = 0.01.*0.15.*0.9.*r.*dip_pos./sqrt(sum(dip_pos.^2,2));
%dip_mom = [0 1 0; 0 0 -1; -1 0 0 ];
dip_pos = pos(segprob.brain(:),:);
dip_pos(sqrt(sum(dip_pos.^2,2))<=1, :) = [];

if size(dip_pos,1) > 100
  dip_pos = dip_pos(randperm(size(dip_pos,1), 100), :);
end

dip_pos = dip_pos./100;
dip_mom = 2.*rand(size(dip_pos)) - 1;

if doplot
  figure
  ft_plot_ortho(segindx.tissue, 'transform', segindx.transform, 'location', center, 'style', 'intersect');
  hold on
  quiver3(dip_pos(:,1),dip_pos(:,2),dip_pos(:,3),dip_mom(:,1),dip_mom(:,2),dip_mom(:,3),'bo')
  view(-200, 15)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. compute the leadfield with singlesphere model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare headmodel
cfg               = [];
cfg.method        = 'singlesphere';
cfg.tissue        = 'brain';
vol_ssph          = ft_prepare_headmodel(cfg, segprob);
vol_ssph          = ft_convert_units(vol_ssph, 'm');

cfg.method        = 'concentricspheres';
cfg.conductivity  = [0.33, 0.43, 0.53];  % vector, conductivity values for tissues: check the order here
cfg.tissue        = {'brain' 'skull' 'scalp'};
vol_csph          = ft_prepare_headmodel(cfg, segprob);
vol_csph          = ft_convert_units(vol_csph, 'm');

%% prepare sourcemodel
cfg                    = [];
cfg.sourcemodel.pos    = dip_pos;
cfg.sourcemodel.inside = ones(size(dip_pos,1),1);
cfg.grad               = sens;
cfg.headmodel          = vol_ssph;
sm_ssph                = ft_prepare_sourcemodel(cfg);

%% prepare leadfield MEG
cfg                 = [];
cfg.sourcemodel     = sm_ssph;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_ssph;
cfg.grad            = sens;
cfg.reducerank      = 3;
out_ssph            = ft_prepare_leadfield(cfg);
lf_meg_ssph         = cell2mat(out_ssph.leadfield);

%% prepare leadfield EEG
cfg                 = [];
cfg.sourcemodel     = sm_ssph;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_csph;
cfg.elec            = elec;
cfg.reducerank      = 3;
out_ssph            = ft_prepare_leadfield(cfg);
lf_eeg_csph         = cell2mat(out_ssph.leadfield);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. compute the leadfield with the duneuro mex-file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare headmodel
cfg              = [];
cfg.method       = 'duneuro';
cfg.conductivity = [0.33, 0.43, 0.53];  % vector, conductivity values for tissues: check the order here
vol_duneuro_hex  = ft_prepare_headmodel(cfg, mesh_vol_hex);
vol_duneuro_tet  = ft_prepare_headmodel(cfg, mesh_vol_tet);

%% prepare sourcemodel
cfg                    = [];
cfg.sourcemodel.pos    = dip_pos;
cfg.sourcemodel.inside = ones(size(dip_pos,1),1);
cfg.grad               = sens;

% hex
cfg.headmodel          = vol_duneuro_hex;
sm_duneuro_hex         = ft_prepare_sourcemodel(cfg);

% tet
cfg.headmodel          = vol_duneuro_tet;
sm_duneuro_tet         = ft_prepare_sourcemodel(cfg);

%% prepare leadfield MEG
cfg                 = [];
cfg.grad            = sens;
cfg.reducerank      = 3;

% hex
cfg.sourcemodel     = sm_duneuro_hex;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_duneuro_hex;
out_hex             = ft_prepare_leadfield(cfg);
lf_meg_hex          = cell2mat(out_hex.leadfield);

% tet
cfg.sourcemodel     = sm_duneuro_tet;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_duneuro_tet;
out_tet             = ft_prepare_leadfield(cfg);
lf_meg_tet          = cell2mat(out_tet.leadfield);

%% prepare leadfield EEG
cfg                 = [];
cfg.elec            = elec;
cfg.reducerank      = 3;

% hex
cfg.sourcemodel     = sm_duneuro_hex;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_duneuro_hex;
out_hex             = ft_prepare_leadfield(cfg);
lf_eeg_hex          = cell2mat(out_hex.leadfield);

% tet
cfg.sourcemodel     = sm_duneuro_tet;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_duneuro_tet;
out_tet             = ft_prepare_leadfield(cfg);
lf_eeg_tet          = cell2mat(out_tet.leadfield);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. compute the leadfield with the brainstorm app
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p,f,e] = fileparts(which('duneuro_matlab'));

% this is needed at the DCCN HPC
setenv('LD_LIBRARY_PATH',[fullfile(p, 'lib') ':' getenv('LD_LIBRARY_PATH')]);

cfg                     = [];
cfg.duneuro.application = duneuro_installbrainstormapp;
cfg.method              = 'duneuro';
cfg.conductivity        = [0.33, 0.43, 0.53];  % vector, conductivity values for tissues: check the order here
vol_duneuro2_hex        = ft_prepare_headmodel(cfg, mesh_vol_hex);
vol_duneuro2_tet        = ft_prepare_headmodel(cfg, mesh_vol_tet);

cfg                 = [];
cfg.grad            = sens;
cfg.reducerank      = 3;

% hexmeg
cfg.sourcemodel     = sm_duneuro_hex;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_duneuro2_hex;
out_hex2            = ft_prepare_leadfield(cfg);
lf_meg_hex2         = cell2mat(out_hex2.leadfield);

% tetmeg
cfg.sourcemodel     = sm_duneuro_tet;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_duneuro2_tet;
out_tet2            = ft_prepare_leadfield(cfg);
lf_meg_tet2         = cell2mat(out_tet2.leadfield);

cfg = rmfield(cfg, 'grad');
cfg.elec = elec;

% hexeeg
cfg.sourcemodel     = sm_duneuro_hex;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_duneuro2_hex;
out_hex2            = ft_prepare_leadfield(cfg);
lf_eeg_hex2         = cell2mat(out_hex2.leadfield);

% teteeg
cfg.sourcemodel     = sm_duneuro_tet;
cfg.sourcemodel.mom = dip_mom';
cfg.headmodel       = vol_duneuro2_tet;
out_tet2            = ft_prepare_leadfield(cfg);
lf_eeg_tet2         = cell2mat(out_tet2.leadfield);

%% compare leadfields
err1 = 100*(abs(lf_meg_hex-lf_meg_ssph)./sqrt(sum(lf_meg_ssph.^2,1)));
err2 = 100*(abs(lf_meg_tet-lf_meg_ssph)./sqrt(sum(lf_meg_ssph.^2,1)));
err3 = 100*(abs(lf_meg_hex2-lf_meg_ssph)./sqrt(sum(lf_meg_ssph.^2,1)));
err4 = 100*(abs(lf_meg_tet2-lf_meg_ssph)./sqrt(sum(lf_meg_ssph.^2,1)));

err = log10([err1(:) err2(:) err3(:) err4(:)]);
hdat = histc(err, -2:0.2:2);
figure;plot(-2:0.2:2, hdat); title('log10 % of meg mismatch with spherical model');

% assert(max(err1(:))<10);
% assert(max(err2(:))<10);
% assert(max(err3(:))<10);
% assert(max(err4(:))<10);

figure; hold on;
plot(lf_meg_hex,'r');
plot(lf_meg_tet,'b');
plot(lf_meg_ssph,'g');
plot(lf_meg_hex2,'k');
plot(lf_meg_tet2,'m');

err1 = 100*(abs(lf_eeg_hex-lf_eeg_csph)./sqrt(sum(lf_eeg_csph.^2,1)));
err2 = 100*(abs(lf_eeg_tet-lf_eeg_csph)./sqrt(sum(lf_eeg_csph.^2,1)));
err3 = 100*(abs(lf_eeg_hex2-lf_eeg_csph)./sqrt(sum(lf_eeg_csph.^2,1)));
err4 = 100*(abs(lf_eeg_tet2-lf_eeg_csph)./sqrt(sum(lf_eeg_csph.^2,1)));

err = log10([err1(:) err2(:) err3(:) err4(:)]);
hdat = histc(err, -2:0.2:2);
figure;plot(-2:0.2:2, hdat); title('log10 % of eeg mismatch with spherical model');

% assert(max(err1(:))<10);
% assert(max(err2(:))<10);
% assert(max(err3(:))<10);
% assert(max(err4(:))<10);





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
%cfg.shrinkSourceSpace = 0;

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


