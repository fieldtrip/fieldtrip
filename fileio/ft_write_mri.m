function [V] = ft_write_mri(filename, dat, varargin)

% FT_WRITE_MRI exports volumetric data such as anatomical and functional MRI
% to a file.
%
% Use as
%   ft_write_mri(filename, dat, ...)
% where the input argument dat represents the 3D array with the values.
%
% The 3-D array with the values can be further described with
%   'transform'    = 4x4 homogenous transformation matrix, specifying the transformation from voxel coordinates to head or world coordinates
%   'unit'         = string, desired geometrical units for the data, for example 'mm'
%   'coordsys'     = string, desired coordinate system for the data
%
% Additional options should be specified in key-value pairs and can be
%   'dataformat'   = string, see below
%   'spmversion'   = string, version of SPM to be used (default = 'spm12')
%   'scl_slope'    = slope parameter for nifti files
%   'scl_inter'    = intersect parameter for nifti files
%   'vmpversion'   = 1 or 2, version of the vmp format to use (default = 2)
%
% The specified filename can already contain the filename extention. If not present,
% it will be added automatically.
%
% The supported dataformats are
%   'analyze'     outdated format and not recommended
%   'mgz'         FreeSurfer specific format
%   'mgh'         FreeSurfer specific format
%   'nifti'     	uses FreeSurfer code
%   'nifti2'      uses FreeSurfer code
%   'nifti_gz'    uses FreeSurfer code
%   'nifti_spm'   uses SPM
%   'vista'       SIMBIO specific format
%   'vmr'         Brainvoyager specific format
%   'vmp'         Brainvoyager specific format
%
% See also FT_READ_MRI, FT_WRITE_DATA, FT_WRITE_HEADSHAPE

% Copyright (C) 2011-2021, Jan-Mathijs Schoffelen
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
transform     = ft_getopt(varargin, 'transform');   % this complements the input data when specified as a 3D array
unit          = ft_getopt(varargin, 'unit');        % this complements the input data when specified as a 3D array
coordsys      = ft_getopt(varargin, 'coordsys');    % this complements the input data when specified as a 3D array
spmversion    = ft_getopt(varargin, 'spmversion');
dataformat    = ft_getopt(varargin, 'dataformat');
scl_slope     = ft_getopt(varargin, 'scl_slope', 1);
scl_inter     = ft_getopt(varargin, 'scl_inter', 0);

% use the version that is on the path, or default to spm12
if ~ft_hastoolbox('spm') && isempty(spmversion)
  spmversion = 'spm12';
end

%% ensure that the input data is consistent
if isnumeric(dat) || islogical(dat)
  % convert the data to a structure according to FT_DATATYPE_VOLUME and FT_READ_MRI
  mri.anatomy = dat;
  mri.dim     = size(dat);
  clear dat % to avoid confusion
  if ~isempty(transform)
    mri.transform = transform;
  else
    mri.coordsys = eye(4);
  end
  if ~isempty(unit)
    mri.unit = unit;
  else
    mri = ft_determine_units(mri);
  end
  if ~isempty(coordsys)
    mri.coordsys = coordsys;
  else
    mri.coordsys = 'unknown';
  end
elseif ft_datatype(dat, 'volume')
  % assume that the data is specified according to FT_DATATYPE_VOLUME and FT_READ_MRI
  mri = dat;
  clear dat % to avoid confusion
  if ~isempty(transform) && ~isequal(transform, mri.transform)
    ft_error('the transformation matrices are not consistent');
  end
  if ~isempty(unit)
    % convert it to the specified units, this will also estimate the units if required
    mri = ft_convert_units(mri, unit);
  end
  if ~isempty(coordsys)
    % convert it to the specified coordinate system, this will also interactively determine the coordinate system if required
    mri = ft_convert_coordsys(mri, coordsys);
  elseif ~isfield(mri, 'coordsys')
    mri.coordsys = 'unknown';
  end
else
  ft_error('unsupported input data');
end

% decompose the input data into separate fields
dat         = mri.anatomy;
transform   = mri.transform;
clear mri % to avoid confusion

%% determine the details of the output format
if isempty(dataformat)
  % only do the autodetection if the format was not specified
  dataformat = ft_filetype(filename);
end

% determine the filename extension
[p, f, x] = fileparts(filename);
if isempty(x)
  switch  dataformat
    case {'vmr', 'vmp', 'mgz', 'mgh'}
      x = ['.' dataformat];
    case {'analyze_old', 'analyze_img', 'analyze_hdr', 'analyze', 'nifti_img'}
      x = '.img';
    case {'freesurfer_mgz'}
      x = '.mgz';
    case {'nifti_spm', 'nifti',  'nifti2'}
      x = '.nii';
    case {'nifti_gz'}
      x = '.nii.gz';
  end
  % ensure that the filename has an extension
  filename = fullfile(p, [f x]);
end

if nargout>0
  % this will be updated by the SPM code and otherwise not used
  V = [];
end

% ensure that the directory exists if we want to write to a file
if ~ismember(dataformat, {'empty', 'fcdc_global', 'fcdc_buffer', 'fcdc_mysql'})
  isdir_or_mkdir(fileparts(filename));
end

%% write the data
switch dataformat
  
  case {'vmr' 'vmp'}
    %% write to brainvoyager format
    vmpversion = ft_getopt(varargin, 'vmpversion', 2);
    write_brainvoyager(filename, dat, dataformat, vmpversion);
    
  case {'analyze_old'}
    %% write to Analyze format, using old code and functions from Darren Webbers toolbox
    avw = avw_hdr_make;
    
    % specify the image data and dimensions
    avw.hdr.dime.dim(2:4) = size(dat);
    avw.img = dat;
    
    % orientation 0 means transverse unflipped (axial, radiological)
    % X direction first,  progressing from patient right to left,
    % Y direction second, progressing from patient posterior to anterior,
    % Z direction third,  progressing from patient inferior to superior.
    avw.hdr.hist.orient = 0;
    
    % specify voxel size
    avw.hdr.dime.pixdim(2:4) = [1 1 1];
    % FIXME, this currently does not work due to all flipping and permuting
    % resx = x(2)-x(1);
    % resy = y(2)-y(1);
    % resz = z(2)-z(1);
    % avw.hdr.dime.pixdim(2:4) = [resy resx resz];
    
    % specify the data type
    switch class(dat)
      case 'logical'
        avw.hdr.dime.datatype = 1;
        avw.hdr.dime.bitpix   = 1;
      case 'uint8'
        avw.hdr.dime.datatype = 2;
        avw.hdr.dime.bitpix   = 8;
      case 'int16'
        avw.hdr.dime.datatype = 4;
        avw.hdr.dime.bitpix   = 16;
      case 'int32'
        avw.hdr.dime.datatype = 8;
        avw.hdr.dime.bitpix   = 32;
      case 'single'
        avw.hdr.dime.datatype = 16;
        avw.hdr.dime.bitpix   = 32;
      case 'double'
        avw.hdr.dime.datatype = 64;
        avw.hdr.dime.bitpix   = 64;
      otherwise
        ft_error('unsupported datatype "%s" to write to analyze_old format', class(dat));
    end
    
    % write the header and image data
    avw_img_write(avw, filename, [], 'ieee-le');
    
  case {'analyze_img' 'analyze_hdr' 'analyze' 'nifti_img' 'nifti_spm'}
    %% analyze or nifti data, using SPM
    V = volumewrite_spm(filename, dat, transform, spmversion, scl_slope, scl_inter);
    
  case {'freesurfer_mgz' 'mgz' 'mgh'}
    %% mgz data, using Freesurfer
    ft_hastoolbox('freesurfer', 1);
    
    % in MATLAB the transformation matrix assumes the voxel indices to be 1-based
    % freesurfer assumes the voxel indices to be 0-based
    transform = vox2ras_1to0(transform);
    save_mgh(dat, filename, transform);
    
  case {'nifti' 'nifti_gz' 'nifti2'}
    %% nifti data, using Freesurfer
    ft_hastoolbox('freesurfer', 1);
    
    datatype = class(dat);
    switch(datatype)
      case 'int8'
        datatype = 'char';
      case 'int16'
        datatype = 'short';
      case 'int32'
        datatype = 'int';
      case 'double'
        datatype = 'double';
      case 'single'
        datatype = 'float';
      case 'uint8'
        datatype = 'uchar';
      otherwise
        ft_error('unsupported datatype "s" to write to nifti format', datatype);
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
    mri.scl_slope = scl_slope;
    mri.scl_inter = scl_inter;
    MRIwrite(mri, filename, datatype);
    
  case {'vista'}
    %% this requires the SIMBIO/Vista toolbox
    ft_hastoolbox('simbio', 1);
    write_vista_vol(size(dat), dat, filename);
    
  otherwise
    ft_error('unsupported format "%s"', dataformat);
end % switch dataformat
