function M = dicom2transform(dcmheader)

% DICOM2TRANSFORM converts the DICOM header parameters into a 4x4 homogenous
% transformation matrix that maps voxel indices to the Patient Coordinate System.
% Note that voxel indices are to be counted starting from 1 (MATLAB and Fortran
% convention, not C/C++ and Python convention). This implementation is known to
% result in a different transformation than FreeSurfer, but corresponds to Horos.
%
% Use as
%   M = dicom2transform(dcmheader)
% where the input argument dcmheader is a structure array with header information for
% each slice. The first structure in the DICOM header array must correspond to slice
% 1 and the last one to slice N.
%
% The header structure for each of the slices must contain
%   dcmheader(i).ImagePositionPatient
%   dcmheader(i).ImageOrientationPatient
%
% The output argument M is a 4x4 homogenous transformation matrix that maps voxel
% indices onto PCS world coordinates in millimeter.
%
% Here are some usefull DICOM references
%   https://doi.org/10.1016/j.jneumeth.2016.03.001
%   https://dicom.innolitics.com/ciods/mr-image/image-plane/00200032
%   https://horosproject.org
%
% See also DCMINFO, LOAD_DICOM_SERIES

% Copyright (C) 2021, Robert Oostenveld
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

rx = dcmheader(1).ImageOrientationPatient(1);
ry = dcmheader(1).ImageOrientationPatient(2);
rz = dcmheader(1).ImageOrientationPatient(3);
cx = dcmheader(1).ImageOrientationPatient(4);
cy = dcmheader(1).ImageOrientationPatient(5);
cz = dcmheader(1).ImageOrientationPatient(6);

% The first value in PixelSpacing is the row spacing, the second value is the column spacing.
vr =  dcmheader(1).PixelSpacing(1);
vc =  dcmheader(1).PixelSpacing(2);

x1 =  dcmheader(1).ImagePositionPatient(1);
y1 =  dcmheader(1).ImagePositionPatient(2);
z1 =  dcmheader(1).ImagePositionPatient(3);

n  = length(dcmheader);
xn =  dcmheader(n).ImagePositionPatient(1);
yn =  dcmheader(n).ImagePositionPatient(2);
zn =  dcmheader(n).ImagePositionPatient(3);

% see equation 1 in https://doi.org/10.1016/j.jneumeth.2016.03.001
Rdicom = [
  rx*vr cx*vc (xn-x1)/(n-1) x1
  ry*vr cy*vc (yn-y1)/(n-1) y1
  rz*vr cz*vc (zn-z1)/(n-1) z1
  0     0      0            1
  ];

% convert from 0-offset to 1-offset voxel locations
Roffset = [
  1 0 0 -1
  0 1 0 -1
  0 0 1 -1
  0 0 0  1
  ];

M = Rdicom*Roffset;
