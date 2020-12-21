function [V] = ft_write_mri(filename, dat, varargin)

% FT_WRITE_MRI exports volumetric data such as anatomical and functional
% MRI to a file.
%
% Use as
%   ft_write_mri(filename, dat, ...)
% where the input argument dat represents the 3-D array with the values.
%
% Additional options should be specified in key-value pairs and can be
%   'dataformat'   = string, see below
%   'transform'    = 4x4 homogenous transformation matrix, specifying the transformation from voxel coordinates to head coordinates
%   'unit'         = string, desired units for the image data on disk, for example 'mm'
%   'spmversion'   = version of SPM to be used, in case data needs to be written in analyze format
%   'scl_slope'    = slope parameter for nifti files
%   'scl_inter'    = intersect parameter for nifti files
%
% The specified filename can already contain the filename extention, but that is not
% required since it will be added automatically.
%
% The supported dataformats are
%   'analyze'
%   'nifti'
%   'vista'
%   'mgz'   (freesurfer)
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
transform     = ft_getopt(varargin, 'transform',  eye(4));
spmversion    = ft_getopt(varargin, 'spmversion', 'spm12');
dataformat    = ft_getopt(varargin, 'dataformat'); % FIXME this is inconsistent with ft_read_mri, which uses 'format'
unit          = ft_getopt(varargin, 'unit');
scl_slope     = ft_getopt(varargin, 'scl_slope', 1);
scl_inter     = ft_getopt(varargin, 'scl_inter', 0);

% convert the input to the desired units
if ~isempty(unit)
  % organize the input data as a FieldTrip structure and estimate its units
  tmp.anatomy   = dat;
  tmp.dim       = size(dat);
  tmp.transform = transform;
  % convert  the input data to the desired units
  tmp = ft_convert_units(tmp, unit);
  % the transformation matrix is the only thing that would have changed
  transform = tmp.transform;
end

if isstruct(dat) && isfield(dat, 'anatomy') && isequal(transform, eye(4))
  % this is an anatomical MRI as returned by FT_READ_MRI
  transform = dat.transform;
  dat       = dat.anatomy;
end

if isempty(dataformat)
  % only do the autodetection if the format was not specified
  dataformat = ft_filetype(filename);
end

if nargout>0
  % start with an empty output argument, it will only be returned by the SPM formats
  V = [];
end

% ensure that the directory exists if we want to write to a file
if ~ismember(dataformat, {'empty', 'fcdc_global', 'fcdc_buffer', 'fcdc_mysql'})
  isdir_or_mkdir(fileparts(filename));
end

switch dataformat
  
  case {'vmr' 'vmp'}
    % brainvoyager file formats
    vmpversion = ft_getopt(varargin, 'vmpversion', 2);
    write_brainvoyager(filename, dat, dataformat, vmpversion);
    
  case {'analyze_old'}
    % analyze format, using old code
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in Analyze format, using some functions from Darren Webbers toolbox
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        ft_error('unsupported datatype %s to write to analyze_old format', class(dat));
    end

    % write the header and image data
    avw_img_write(avw, filename, [], 'ieee-le');
    
  case {'analyze_img' 'analyze_hdr' 'analyze' 'nifti_img' 'nifti_spm'}
    % analyze or nifti data, using SPM
    V = volumewrite_spm(filename, dat, transform, spmversion, scl_slope, scl_inter);
    
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
        ft_error('unsupported datatype to write to a nifti file');
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
    % this requires the SIMBIO/Vista toolbox
    ft_hastoolbox('simbio', 1);
    write_vista_vol(size(dat), dat, filename);
    
  otherwise
    ft_error('unsupported format "%s"', dataformat);
end % switch dataformat
