% Copyright (C) 2004 Yokogawa Electric Corporation, All Rights Reserved.
%
% MRI data loading program
% 
% USAGE:
%       [ header_info, img_data ] = GetMeg160MriFileM( fid );
%
% ARGUMENT:
%   fid             : file ID
%
% RETURN VALUES:
%   header_info : Structure which has fields about .mri header information.
%   img_data    : Structure 
%
function [ header_info, view_count, view_list ] = GetMeg160MriFileM( filename );
