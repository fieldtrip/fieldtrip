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
strPath_Data = '/home/common/matlab/fieldtrip/data/test/bug2336/';


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


% Read electrode config, surface points & fiducials, FIXME this function is missing
% sfh = readBESAsfh(strPath__ElectrodeConfig);


% Read scalp, skull, csf & brain surface
[pnt_scalp, tri_scalp, srf_scalp] = ft_read_bv_srf(strPath__ScalpSurface);
[pnt_skull, tri_skull, srf_skull] = ft_read_bv_srf(strPath__SkullSurface);
[pnt_csf, tri_csf, srf_csf]       = ft_read_bv_srf(strPath__CSFSurface);
[pnt_brain, tri_brain, srf_brain] = ft_read_bv_srf(strPath__BrainSurface);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The second part of the test script pertains to the high-level fieldtrip implementation

cfg = [];
cfg.method = 'besa';
cfg.elecfile = strPath__ElectrodeConfig;
cfg.hdmfile = strPath_Leadfield;
cfg.outputfile = tempname;
cfg.smooth = 'no';
vol = ft_prepare_headmodel(cfg);
