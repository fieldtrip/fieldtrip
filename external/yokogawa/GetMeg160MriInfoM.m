% Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
%
% usage:
%   mri_info = GetMeg160MriInfoM( fid )
%
% arguments:
%   fid             : file ID
%
% return values:
%   mri_info        : Structure of MRI information.
%                       .type_name      : MRI type name ('NoMriFile', 'NormalMriFile', 'VirtualMriFile')
%                       .file_name      : MRI file name
%                       .model_name     : conductor model type ('NoModel', 'SphericalModel')
%                       .cx, .cy, .cz   : spherical center position for 'SphericalModel'
%                       .r              : spherical radius for 'SphericalModel'
function    mri_info  =   GetMeg160MriInfoM( fid )