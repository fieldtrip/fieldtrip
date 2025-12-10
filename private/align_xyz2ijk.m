function [volume, permutevec, flipvec] = align_xyz2ijk(volume)

% ALIGN_XYZ2IJK updates the transform and coordsys fields such that the axes of the
% resulting head coordinate system are aligned with the voxel indices. The intention
% is to create a volume structure that can be plotted in native voxel coordinates.
%
% See also ALIGN_IJK2XYZ, VOLUMEPERMUTE, VOLUMEFLIP

% Copyright (C) 2023, Jan-Mathijs Schoffelen and Robert Oostenveld
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

% convert 'neuromag' to 'ras', etc.
volume.coordsys = generic(volume.coordsys);

% permute the 3D volume so that the indices and head-coordinate axes correspond approximately
% the determination of the permutation order requires that voxel size differneces are accounted for
voxsiz = sqrt(sum(volume.transform(:,1:3).^2,1));
[dum, dim1] = max(abs(volume.transform(1,1:3))./voxsiz);
[dum, dim2] = max(abs(volume.transform(2,1:3))./voxsiz);
[dum, dim3] = max(abs(volume.transform(3,1:3))./voxsiz);
permutevec  = [dim1 dim2 dim3];

if length(unique(permutevec))<3
  ft_error('could not determine the correspondence between volume and headcoordinate axes');
end

% permute the axes of the transform and coordsys
volume.transform(permutevec,:) = volume.transform(1:3,:);
volume.coordsys(permutevec) = volume.coordsys;

% flip the volume along each direction, so that the diagonal of the transformation matrix is positive
flipvec = diag(volume.transform(1:3,1:3))<0;
if flipvec(1)
  volume.transform(1,:) = -1 * volume.transform(1,:);
  volume.coordsys(1) = flipletter(volume.coordsys(1));
end
if flipvec(2)
  volume.transform(2,:) = -1 * volume.transform(2,:);
  volume.coordsys(2) = flipletter(volume.coordsys(2));
end
if flipvec(3)
  volume.transform(3,:) = -1 * volume.transform(3,:);
  volume.coordsys(3) = flipletter(volume.coordsys(3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function letter = flipletter(letter)
switch letter
  case 'a'
    letter = 'p';
  case 'p'
    letter = 'a';
  case 'l'
    letter = 'r';
  case 'r'
    letter = 'l';
  case 'i'
    letter = 's';
  case 's'
    letter = 'i';
  otherwise
    ft_error('incorrect letter')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coordsys = generic(coordsys)
mapping = {
  'ctf',       'als'
  'bti',       'als'
  '4d',        'als'
  'yokogawa',  'als'
  'eeglab',    'als'
  'neuromag',  'ras'
  'itab',      'ras'
  'acpc',      'ras'
  'spm',       'ras'
  'mni',       'ras'
  'fsaverage', 'ras'
  'tal',       'ras'
  'scanras',   'ras'
  'scanlps',   'lps'
  'dicom',     'lps'
  'paxinos',   'rsp'
  };

sel = find(strcmp(mapping(:,1), coordsys));
if length(sel)==1
  coordsys = mapping{sel,2};
end

if ~all(ismember(coordsys, 'lrapis'))
  ft_error('cannot convert "%s" to a generic coordinate system label', coordsys);
end
