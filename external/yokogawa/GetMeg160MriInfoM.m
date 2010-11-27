%  Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
% 
%  usage:
%    mri_info = GetMeg160MriInfoM( fid )
% 
%  arguments:
%    fid             : file ID
% 
%  return values:
%    mri_info        : Structure of MRI information.
%                        .type_name      : MRI type name ('NoMriFile', 'NormalMriFile', 'VirtualMriFile')
%                        .file_name      : MRI file name
%                        .model_name     : conductor model type ('NoModel', 'SphericalModel')
%                        .cx, .cy, .cz   : spherical center position for 'SphericalModel' [m]
%                        .r              : spherical radius for 'SphericalModel' [m]
%  
%  confirmation of revision:
%   GetMeg160MriInfoM( Inf ) will show and return revision of this function.
%
