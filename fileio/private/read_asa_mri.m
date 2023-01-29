function [mri, seg, hdr] = read_asa_mri(fn)

% READ_ASA_MRI reads an ASA format MRI file
%
% Use as
%   [mri, seg, hdr] = read_asa_mri(filename)
%
% The raw image data is returned, together with the position of the
% external head markers in raw image coordinates.
%
% In the ASA default PAN (pre-auricular/nasion) coordinate system
%   PointOnPositiveYAxis -> LPA
%   PointOnNegativeYAxis -> RPA
%   PointOnPositiveXAxis -> nasion

% Copyright (C) 2002, Robert Oostenveld
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

hdr.Nrows     = read_ini(fn, 'NumberRows=', '%d');
hdr.Ncolumns  = read_ini(fn, 'NumberColumns=', '%d');
hdr.Nslices   = read_ini(fn, 'NumberSlices=', '%d');
hdr.rows      = read_ini(fn, 'Rows=', '%s');
hdr.columns   = read_ini(fn, 'Columns=', '%s');
hdr.slices    = read_ini(fn, 'Slices=', '%s');
hdr.distance  = read_ini(fn, 'Distance=', '%f', 3);
hdr.mrifile   = read_ini(fn, 'Matrix', '%s');
hdr.segfile   = read_ini(fn, 'Segmentation', '%s');
hdr.posy      = read_ini(fn, 'PointOnPositiveYAxis', '%f', 3);
hdr.negy      = read_ini(fn, 'PointOnNegativeYAxis', '%f', 3);
hdr.posx      = read_ini(fn, 'PointOnPositiveXAxis', '%f', 3);
hdr.voxposy   = read_ini(fn, 'VoxelOnPositiveYAxis', '%f', 3);
hdr.voxnegy   = read_ini(fn, 'VoxelOnNegativeYAxis', '%f', 3);
hdr.voxposx   = read_ini(fn, 'VoxelOnPositiveXAxis', '%f', 3);
hdr.segfile   = read_ini(fn, 'Segmentation', '%s', 1);

dim = [hdr.Ncolumns hdr.Nrows hdr.Nslices];
mri = [];
seg = [];

% the data files are at the same location as the mri file, locate path
[path, name, ext] = fileparts(fn);

% this temporary needs 8x as much storage!!!
if ~isempty(hdr.mrifile)
  mrifile = fullfile(path, hdr.mrifile);
  mri = zeros(dim);
  fid = fopen_or_error(mrifile, 'rb', 'ieee-le');
  mri(:) = fread(fid, prod(dim), 'uint8');
  mri = uint8(mri);
  fclose(fid);
end

% this temporary needs 8x as much storage!!!
if ~isempty(hdr.segfile)
  segfile = fullfile(path, hdr.segfile);
  seg = zeros(dim);
  fid = fopen_or_error(segfile, 'rb', 'ieee-le');
  seg(:) = fread(fid, prod(dim), 'uint8');
  seg = uint8(seg);
  fclose(fid);
end

% flip the orientation of the MRI data, the result should be
% 'coronal(occipital-frontal)'
% 'horizontal(inferior-superior)'
% 'sagittal(right-left)'

if strcmp(hdr.columns, 'coronal(frontal-occipital)')
  hdr.columns = 'coronal(occipital-frontal)';
  mri = flipdim(mri, 1);
  seg = flipdim(seg, 1);
elseif strcmp(hdr.columns, 'horizontal(superior-inferior)')
  hdr.columns = 'horizontal(inferior-superior)';
  mri = flipdim(mri, 1);
  seg = flipdim(seg, 1);
elseif strcmp(hdr.columns, 'sagittal(left-right)')
  hdr.columns = 'sagittal(right-left)';
  mri = flipdim(mri, 1);
  seg = flipdim(seg, 1);
end

if strcmp(hdr.rows, 'coronal(frontal-occipital)')
  hdr.rows = 'coronal(occipital-frontal)';
  mri = flipdim(mri, 2);
  seg = flipdim(seg, 2);
elseif strcmp(hdr.rows, 'horizontal(superior-inferior)')
  hdr.rows = 'horizontal(inferior-superior)';
  mri = flipdim(mri, 2);
  seg = flipdim(seg, 2);
elseif strcmp(hdr.rows, 'sagittal(left-right)')
  hdr.rows = 'sagittal(right-left)';
  mri = flipdim(mri, 2);
  seg = flipdim(seg, 2);
end

if strcmp(hdr.slices, 'coronal(frontal-occipital)')
  hdr.slices = 'coronal(occipital-frontal)';
  mri = flipdim(mri, 3);
  seg = flipdim(seg, 3);
elseif strcmp(hdr.slices, 'horizontal(superior-inferior)')
  hdr.slices = 'horizontal(inferior-superior)';
  mri = flipdim(mri, 3);
  seg = flipdim(seg, 3);
elseif strcmp(hdr.slices, 'sagittal(left-right)')
  hdr.slices = 'sagittal(right-left)';
  mri = flipdim(mri, 3);
  seg = flipdim(seg, 3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% swap the orientations of the MRI data, the result should be fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st dimension corresponds to columns, which should be 'coronal(occipital-frontal)'
% 2st dimension corresponds to rows,    which should be 'sagittal(right-left)'
% 3rd dimension corresponds to slices,  which should be 'horizontal(inferior-superior)'

if     strcmp(hdr.columns, 'coronal(occipital-frontal)')
  orientation(1) = 1;
elseif strcmp(hdr.columns, 'sagittal(right-left)')
  orientation(1) = 2;
elseif strcmp(hdr.columns, 'horizontal(inferior-superior)')
  orientation(1) = 3;
end

if     strcmp(hdr.rows, 'coronal(occipital-frontal)')
  orientation(2) = 1;
elseif strcmp(hdr.rows, 'sagittal(right-left)')
  orientation(2) = 2;
elseif strcmp(hdr.rows, 'horizontal(inferior-superior)')
  orientation(2) = 3;
end

if     strcmp(hdr.slices, 'coronal(occipital-frontal)')
  orientation(3) = 1;
elseif strcmp(hdr.slices, 'sagittal(right-left)')
  orientation(3) = 2;
elseif strcmp(hdr.slices, 'horizontal(inferior-superior)')
  orientation(3) = 3;
end

mri = ipermute(mri, orientation);
seg = ipermute(seg, orientation);
hdr.rows    = 'sagittal(right-left)';
hdr.columns = 'coronal(occipital-frontal)';
hdr.slices  = 'horizontal(inferior-superior)';

% recompute the dimensions after all the swapping
hdr.Nrows    = size(mri, 1);
hdr.Ncolumns = size(mri, 2);
hdr.Nslices  = size(mri, 3);
dim = size(mri);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if possible, create the accompanying homogenous coordinate transformation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In case of PointOn..., ASA counts voxels from the center of the MRI
% and in case of VoxelOn..., ASA counts voxels from the corner of the MRI
% In both cases, ASA starts counting at [0 0 0], which is C convention
% whereas I want to count from the 1st voxel and number that with [1 1 1]
if ~isempty(hdr.posx) && ~isempty(hdr.negy) && ~isempty(hdr.posy)
  offset = (dim + [1 1 1])/2;
  hdr.fiducial.mri.nas = hdr.posx + offset;
  hdr.fiducial.mri.lpa = hdr.posy + offset;
  hdr.fiducial.mri.rpa = hdr.negy + offset;
else
  offset = [1 1 1];
  hdr.fiducial.mri.nas = hdr.voxposx + offset;
  hdr.fiducial.mri.lpa = hdr.voxposy + offset;
  hdr.fiducial.mri.rpa = hdr.voxnegy + offset;
end

% use the headcoordinates function (roboos/misc) to compute the transformaton matrix
hdr.transformMRI2Head = ft_headcoordinates(hdr.fiducial.mri.nas, hdr.fiducial.mri.lpa, hdr.fiducial.mri.rpa, 'asa');
hdr.transformHead2MRI = inv(hdr.transformMRI2Head);

% compute the fiducials in head coordinates
hdr.fiducial.head.nas = ft_warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.nas, 'homogenous');
hdr.fiducial.head.lpa = ft_warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.lpa, 'homogenous');
hdr.fiducial.head.rpa = ft_warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.rpa, 'homogenous');


