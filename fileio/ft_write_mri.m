function [V] = ft_write_mri(filename, dat, varargin)

% FT_WRITE_MRI exports volumetric data such as anatomical and functional
% MRI to a file.
%
% Use as
%   V = ft_write_mri(filename, dat, ...)
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
%   vista
%   mgz   (freesurfer)
%
% See also FT_READ_MRI, FT_WRITE_DATA, FT_WRITE_HEADSHAPE

% Copyright (C) 2011-2012, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
transform     = ft_getopt(varargin, 'transform', eye(4));
spmversion    = ft_getopt(varargin, 'spmversion', 'SPM8');
dataformat    = ft_getopt(varargin, 'dataformat'); % FIXME this is inconsistent with ft_read_mri, which uses 'format'

if isempty(dataformat)
  % only do the autodetection if the format was not specified
  dataformat = ft_filetype(filename);
end

if nargout>0
  % start with an empty output argument, it will only be returned by the SPM formats
  V = [];
end

switch dataformat
  case {'analyze_img' 'analyze_hdr' 'analyze' 'nifti_spm'}
    % analyze data, using SPM
    V = volumewrite_spm(filename, dat, transform, spmversion);
    
  case {'freesurfer_mgz' 'mgz' 'mgh'}
    % mgz data, using Freesurfer
    ft_hastoolbox('freesurfer', 1);
    
    % in MATLAB the transformation matrix assumes the voxel indices to be 1-based
    % freesurfer assumes the voxel indices to be 0-based
    transform = vox2ras_1to0(transform);
    save_mgh(dat, filename, transform);
    
  case {'nifti'}
    % nifti data, using Freesurfer
    ft_hastoolbox('freesurfer', 1);
    
    datatype = class(dat);
    switch(datatype)
      case 'int8'
        datatype = 'uchar';
      case 'int16'
        datatype = 'short';
      case 'int32'
        datatype = 'int';
      case 'double'
        datatype = 'double';
      case 'single'
        datatype = 'float';
      case 'logical'
        datatype = 'uchar';
      otherwise
        error('unsupported datatype to write to Nifti');
    end
    
    ndims = numel(size(dat));
    if ndims==3
      dat = ipermute(dat, [2 1 3]);
      % FIXME although this is probably correct
      % see the help of MRIread, anecdotally columns and rows seem to need a swap in
      % order to match the transform matrix (alternatively a row switch of the latter
      % can be done) to keep the writing consistent with the reading
    elseif ndims==4
      dat = ipermute(dat, [2 1 3 4]);
    end
    
    mri          = [];
    mri.vol      = dat;
    mri.vox2ras0 = vox2ras_1to0(transform);
    mri.volres   = sqrt(sum(transform(:,1:3).^2));
    MRIwrite(mri, filename, datatype);
    
  case {'vista'}
    % this requires the SIMBIO/Vista toolbox
    ft_hastoolbox('simbio', 1);
    write_vista_vol(size(dat), dat, filename);
    
  otherwise
    error('unsupported format "%s"', dataformat);
end % switch dataformat
