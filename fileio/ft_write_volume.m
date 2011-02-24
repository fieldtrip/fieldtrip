function [V] = ft_write_volume(filename, dat, varargin)

% FT_WRITE_VOLUME exports volumetric data to a file.
%
% Use as
%   V = ft_write_volume(filename, dat, ...)
%
% The specified filename can already contain the filename extention,
% but that is not required since it will be added automatically.
%
% Additional options should be specified in key-value pairs and can be
%   'spmversion'     spmversion to be used (in case data needs to be
%                      written in analyze format
%   'dataformat'     string, see below
%   'transform'      transformation matrix, specifying the transformation
%                      from voxel coordinates to head coordinates
%
% The supported dataformats are
%   analyze
%   nifti
%   mgz (freesurfer)
%
% See also FT_WRITE_DATA 

% Copyright (C) 2011 Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% get the options
dataformat    = keyval('dataformat',    varargin); if isempty(dataformat), dataformat = ft_filetype(filename); end
transform     = keyval('transform',     varargin); if isempty(transform),  transform  = eye(4);                end
spmversion    = keyval('spmversion',    varargin);

if isempty(spmversion), spmversion = 'SPM8'; end

if strcmp(dataformat, 'nifti') && strcmp(spmversion, 'SPM2') 
  error('nifti can only be written by spm versions newer than spm2');
end

switch dataformat
   
  case {'analyze_img' 'analyze_hdr' 'analyze'}

    %analyze data, using SPM
    V = volumewrite_spm(filename, dat, transform, spmversion);

  case {'freesurfer_mgz' 'mgz'}
    % mgz-volume using freesurfer
    ft_hastoolbox('freesurfer', 1);

    % in matlab the transformation matrix assumes the voxel indices to be 1-based
    % freesurfer assumes the voxel indices to be 0-based
    transform = vox2ras_1to0(transform);  
    save_mgh(dat, filename, transform);
    V = [];
 
  case {'nifti'}
    %nifti data, using SPM
    V = volumewrite_spm(filename, dat, transform, spmversion); 
  
  otherwise
    error('unsupported data format');
end % switch dataformat
