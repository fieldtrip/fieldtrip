function test_bug2336

% these settings are needed to execute the test script on out cluster
% WALLTIME 00:30:00
% MEM 4gb

% script to test the integration of FEM output files created by BESA MRI.
%
% Copyright (C) 2013, BESA GmbH
%
% File name: test_bug2336.m
%
% Author: Robert Spangler
% Created: 2013-11-27


% ensure that fieldtrip/extrenal/besa is on the path
ft_hastoolbox('besa', 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The first part of the test script pertains to the low-level reading functions

% Specify path to data folder
strPath_Data = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2336/');

% Filenames of BESA MRI output files
strPath_Leadfield         = sprintf('%sBESA_MRI_Electrode_Config_34_PB_FEM_DATA.lft', strPath_Data);
strPath_SourceSpaceNodes  = sprintf('%sBESA_MRI_Electrode_Config_34_PB_FEM_DATA.loc', strPath_Data);
strPath__ElectrodeConfig  = sprintf('%sElectrode_Config_34.sfh', strPath_Data);
strPath__SkullSurface     = sprintf('%sBESA_MRI_MRI_PB_FEM_SKULL_ACPC.srf', strPath_Data);
strPath__BrainSurface     = sprintf('%sBESA_MRI_PB_FEM_BRAIN_ACPC.srf', strPath_Data);
strPath__CSFSurface       = sprintf('%sBESA_MRI_PB_FEM_CSF_ACPC.srf', strPath_Data);
strPath__ScalpSurface     = sprintf('%sBESA_MRI_PB_FEM_SCALP_ACPC.srf', strPath_Data);


% Read leadfield
[dim lf] = readBESAlft(strPath_Leadfield);
fprintf('Number of channels: %i\n', dim(1));
fprintf('Number of nodes: %i\n', dim(2));
fprintf('Number of orientations: %i\n', dim(3));


% Read source space grid nodes
[ssg IdxNeighbour] = readBESAloc(strPath_SourceSpaceNodes);
fprintf('Number of nodes: %i\n', size(ssg, 1));
fprintf('Number of neighbours/node: %i\n', size(IdxNeighbour, 1));

% Read electrode config &  fiducials
cfg_eeg = readBESAsfh(strPath__ElectrodeConfig, 'eeg');
cfg_fid = readBESAsfh(strPath__ElectrodeConfig, 'fiducial');

% Read scalp, skull, csf & brain surface
cfg_scalp = readBESAsrf(strPath__ScalpSurface);
cfg_skull = readBESAsrf(strPath__SkullSurface);
cfg_csf   = readBESAsrf(strPath__CSFSurface);
cfg_brain = readBESAsrf(strPath__BrainSurface);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The second part of the test script pertains to the high-level FieldTrip implementation

cfg = [];
cfg.method = 'besa';
cfg.elecfile = strPath__ElectrodeConfig;
cfg.hdmfile = strPath_Leadfield;
cfg.outputfile = tempname;
cfg.smooth = 'no';
vol = ft_prepare_headmodel(cfg);

pos = ssg(1,:);
sens1 = ft_read_sens(strPath__ElectrodeConfig);
sens2 = sens1;
sens2.elecpos = sens2.elecpos+randn(size(sens2.elecpos));

[vol1, sens1] = ft_prepare_vol_sens(vol, sens1);
lf1 = ft_compute_leadfield(pos, sens1, vol1);

[vol2, sens2] = ft_prepare_vol_sens(vol, sens2);
lf2 = ft_compute_leadfield(pos, sens2, vol2);

figure
subplot(2,2,1)
ft_plot_topo3d(sens1.elecpos, lf(:,1))
subplot(2,2,2)
ft_plot_topo3d(sens1.elecpos, lf(:,2))
subplot(2,2,3)
ft_plot_topo3d(sens1.elecpos, lf(:,3))


% show 3D stuff
figure
ft_plot_dipole(pos, [0 0 1]', 'diameter', 10, 'length', 20);
ft_plot_topo3d(sens1.elecpos, lf(:,3))

[X, Y, Z] = ndgrid(1:vol.dim(1), 1:vol.dim(2), 1:vol.dim(3));
pos = ft_warp_apply(vol.transform, [X(:) Y(:) Z(:)]);
ft_plot_mesh(pos);



