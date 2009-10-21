function [mri, seg, hdr] = read_asa_mri(fn);

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
% $Log: read_asa_mri.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.6  2008/11/12 17:02:03  roboos
% explicitely specify ieee-le in fopen()
%
% Revision 1.5  2005/11/16 13:46:32  roboos
% added segmentatino to output, changed from warpo3d to warp_apply
%
% Revision 1.4  2004/01/19 14:24:13  roberto
% numerous changes, cannot remember the details
%
% Revision 1.3  2003/08/04 09:26:46  roberto
% added homgenous coordinate transformation matrices to header
% support for VoxelOn... instead of PointOn...
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

hdr.Nrows     = read_asa(fn, 'NumberRows=', '%d');
hdr.Ncolumns  = read_asa(fn, 'NumberColumns=', '%d');
hdr.Nslices   = read_asa(fn, 'NumberSlices=', '%d');
hdr.rows      = read_asa(fn, 'Rows=', '%s');
hdr.columns   = read_asa(fn, 'Columns=', '%s');
hdr.slices    = read_asa(fn, 'Slices=', '%s');
hdr.distance  = read_asa(fn, 'Distance=', '%f', 3);
hdr.mrifile   = read_asa(fn, 'Matrix', '%s');
hdr.segfile   = read_asa(fn, 'Segmentation', '%s');
hdr.posy      = read_asa(fn, 'PointOnPositiveYAxis', '%f', 3);
hdr.negy      = read_asa(fn, 'PointOnNegativeYAxis', '%f', 3);
hdr.posx      = read_asa(fn, 'PointOnPositiveXAxis', '%f', 3);
hdr.voxposy   = read_asa(fn, 'VoxelOnPositiveYAxis', '%f', 3);
hdr.voxnegy   = read_asa(fn, 'VoxelOnNegativeYAxis', '%f', 3);
hdr.voxposx   = read_asa(fn, 'VoxelOnPositiveXAxis', '%f', 3);
hdr.segfile   = read_asa(fn, 'Segmentation', '%s', 1);

dim = [hdr.Ncolumns hdr.Nrows hdr.Nslices];
mri = [];
seg = [];

% the data files are at the same location as the mri file, locate path
[path, name, ext] = fileparts(fn);

% this temporary needs 8x as much storage!!!
if ~isempty(hdr.mrifile)
  mrifile = fullfile(path, hdr.mrifile);
  mri = zeros(dim);
  fid = fopen(mrifile, 'rb', 'ieee-le');
  mri(:) = fread(fid, prod(dim), 'uint8');
  mri = uint8(mri);
  fclose(fid);
end

% this temporary needs 8x as much storage!!!
if ~isempty(hdr.segfile)
  segfile = fullfile(path, hdr.segfile);
  seg = zeros(dim);
  fid = fopen(segfile, 'rb', 'ieee-le');
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

if exist('headcoordinates', 'file')
  % In case of PointOn..., ASA counts voxels from the center of the MRI
  % and in case of VoxelOn..., ASA counts voxels from the corner of the MRI
  % In both cases, ASA starts counting at [0 0 0], which is C convention
  % whereas I want to count from the 1st voxel and number that with [1 1 1]
  if ~isempty(hdr.posx) & ~isempty(hdr.negy) & ~isempty(hdr.posy)
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
  hdr.transformMRI2Head = headcoordinates(hdr.fiducial.mri.nas, hdr.fiducial.mri.lpa, hdr.fiducial.mri.rpa, 1);
  hdr.transformHead2MRI = inv(hdr.transformMRI2Head);

  % compute the fiducials in head coordinates
  hdr.fiducial.head.nas = warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.nas, 'homogenous');
  hdr.fiducial.head.lpa = warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.lpa, 'homogenous');
  hdr.fiducial.head.rpa = warp_apply(hdr.transformMRI2Head, hdr.fiducial.mri.rpa, 'homogenous');
end

