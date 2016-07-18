function [h, coordsys] = ft_headcoordinates(fid1, fid2, fid3, fid4, coordsys)

% FT_HEADCOORDINATES returns the homogeneous coordinate transformation matrix
% that converts the specified fiducials in any coordinate system (e.g. MRI)
% into the rotated and translated headcoordinate system.
%
% Use as
%   [h, coordsys] = ft_headcoordinates(fid1, fid2, fid3, coordsys)
% or
%   [h, coordsys] = ft_headcoordinates(fid1, fid2, fid3, fid4, coordsys)
%
% Depending on the desired coordinate system, the order of the fiducials is
% interpreted as follows
%
%   fid1 = nas
%   fid2 = lpa
%   fid3 = rpa
%   fid4 = extra point (optional)
%
%   fid1 = ac
%   fid2 = pc
%   fid3 = midsagittal
%   fid4 = extra point (optional)
%
%   fid1 = pt1
%   fid2 = pt2
%   fid3 = pt3
%   fid4 = extra point (optional)
%
%   fid1 = bregma
%   fid2 = lambda
%   fid3 = midsagittal
%   fid4 = extra point (optional)
%
% The fourth argument fid4 is optional and can be specified as an an extra point
% which is assumed to have a positive Z-coordinate. It will be used to ensure correct
% orientation of the z-axis (ctf, 4d, yokogawa, itab, neuromag) or X-axis (tal, spm).
% The specification of this extra point may result in the handedness of the
% transformation to be changed, but ensures consistency with the handedness of the
% input coordinate system.
%
% The coordsys input argument is a string that determines how the location of the
% origin and the direction of the axis is to be defined relative to the fiducials
%   according to CTF conventions:             coordsys = 'ctf'
%   according to 4D conventions:              coordsys = '4d' or 'bti'
%   according to YOKOGAWA conventions:        coordsys = 'yokogawa'
%   according to ASA conventions:             coordsys = 'asa'
%   according to NEUROMAG conventions:        coordsys = 'itab'
%   according to ITAB conventions:            coordsys = 'neuromag'
%   according to FTG conventions:             coordsys = 'ftg'
%   according to Talairach conventions:       coordsys = 'tal'
%   according to SPM conventions:             coordsys = 'spm'
%   according to PAXINOS conventions:         coordsys = 'paxinos'
% If coordsys is not specified, it will default to 'ctf'.
%
% The CTF, 4D and YOKOGAWA coordinate systems are defined as follows:
%   the origin is exactly between lpa and rpa
%   the X-axis goes towards nas
%   the Y-axis goes approximately towards lpa, orthogonal to X and in the plane spanned by the fiducials
%   the Z-axis goes approximately towards the vertex, orthogonal to X and Y
%
% The TALAIRACH and SPM coordinate systems are defined as:
%   the origin corresponds with the anterior commissure
%   the Y-axis is along the line from the posterior commissure to the anterior commissure
%   the Z-axis is towards the vertex, in between the hemispheres
%   the X-axis is orthogonal to the midsagittal-plane, positive to the right
%
% The NEUROMAG and ITAB coordinate systems are defined as follows:
%   the X-axis is from the origin towards the RPA point (exactly through)
%   the Y-axis is from the origin towards the nasion (exactly through)
%   the Z-axis is from the origin upwards orthogonal to the XY-plane
%   the origin is the intersection of the line through LPA and RPA and a line orthogonal to L passing through the nasion
%
% The ASA coordinate system is defined as follows:
%   the origin is at the orthogonal intersection of the line from rpa-lpa and the line trough nas
%   the X-axis goes towards nas
%   the Y-axis goes through rpa and lpa
%   the Z-axis goes approximately towards the vertex, orthogonal to X and Y
%
% The FTG coordinate system is defined as:
%   the origin corresponds with pt1
%   the x-axis is along the line from pt1 to pt2
%   the z-axis is orthogonal to the plane spanned by pt1, pt2 and pt3
%
% The PAXINOS coordinate system is defined as:
%   the origin is at bregma
%   the x-axis extends along the Medial-Lateral direction, with positive towards the right
%   the y-axis points from dorsal to ventral, i.e. from inferior to superior
%   the z-axis passes through bregma and lambda and points from cranial to caudal, i.e. from anterior to posterior
%
% See also FT_ELECTRODEREALIGN, FT_VOLUMEREALIGN, FT_INTERACTIVEREALIGN

% Copyright (C) 2003-2014 Robert Oostenveld
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

% figure out the input arguments
if nargin==3
  coordsys = 'ctf'; % default
  fid4     = [];
elseif nargin==4 && ischar(fid4)
  coordsys = fid4;
  fid4     = [];
elseif nargin==4 && isnumeric(fid4)
  coordsys = 'ctf'; % default
elseif nargin==5
  % do nothing
else
  error('incorrect specification of input parameters');
end

if isnumeric(coordsys)
  % these are for backward compatibility, but should preferably not be used any more
  if coordsys==0,
    coordsys = 'ctf';
  elseif coordsys==1,
    coordsys = 'asa';
  elseif coordsys==2,
    coordsys = 'ftg';
  else
    error('if coordsys is numeric, it should assume one of the values 0/1/2');
  end
end

% ensure that they are row vectors
fid1 = fid1(:)';
fid2 = fid2(:)';
fid3 = fid3(:)';
fid4 = fid4(:)';

assert(numel(fid1)==3, 'incorrect specification of fiducial 1');
assert(numel(fid2)==3, 'incorrect specification of fiducial 2');
assert(numel(fid3)==3, 'incorrect specification of fiducial 3');
assert(isempty(fid4) || numel(fid4)==3, 'incorrect specification of fiducial 4');

% ensure that it is lower case
coordsys = lower(coordsys);

% compute the origin and direction of the coordinate axes in MRI coordinates
switch coordsys
  case {'als_ctf' 'ctf' 'bti' '4d' 'yokogawa'}
    % rename the marker points for convenience
    nas = fid1; lpa = fid2; rpa = fid3; extrapoint = fid4; clear fid*
    origin = (lpa+rpa)/2;
    dirx = nas-origin;
    dirz = cross(dirx,lpa-rpa);
    diry = cross(dirz,dirx);
    dirx = dirx/norm(dirx);
    diry = diry/norm(diry);
    dirz = dirz/norm(dirz);
  case {'als_asa' 'asa'}
    % rename the marker points for convenience
    nas = fid1; lpa = fid2; rpa = fid3; extrapoint = fid4; clear fid*
    dirz = cross(nas-rpa, lpa-rpa);
    diry = lpa-rpa;
    dirx = cross(diry,dirz);
    dirx = dirx/norm(dirx);
    diry = diry/norm(diry);
    dirz = dirz/norm(dirz);
    origin = rpa + dot(nas-rpa,diry)*diry;
  case {'ras_itab' 'itab' 'neuromag'}
    % rename the fiducials
    nas = fid1; lpa = fid2; rpa = fid3; extrapoint = fid4; clear fid*
    dirz = cross(rpa-lpa,nas-lpa);
    dirx = rpa-lpa;
    diry = cross(dirz,dirx);
    dirx = dirx/norm(dirx);
    diry = diry/norm(diry);
    dirz = dirz/norm(dirz);
    origin = lpa + dot(nas-lpa,dirx)*dirx;
  case 'ftg'
    % rename the marker points for convenience
    pt1 = fid1; pt2 = fid2; pt3 = fid3; extrapoint = fid4; clear fid*
    origin = pt1;
    dirx = pt2-origin;
    diry = pt3-origin;
    dirz = cross(dirx,diry);
    diry = cross(dirz,dirx);
    dirx = dirx/norm(dirx);
    diry = diry/norm(diry);
    dirz = dirz/norm(dirz);
  case {'ras_tal' 'tal' 'spm'}
    % rename the marker points for convenience
    ac = fid1; pc = fid2; midsagittal = fid3; extrapoint = fid4; clear fid*
    origin = ac;
    diry   = ac-pc;
    dirz   = midsagittal-ac;
    dirx   = cross(diry,dirz);
    dirz   = cross(dirx,diry);
    dirx   = dirx/norm(dirx);
    diry   = diry/norm(diry);
    dirz   = dirz/norm(dirz);
  case 'paxinos'
    bregma = fid1; lambda = fid2; midsagittal = fid3; extrapoint = fid4; clear fid*
    origin = bregma;           % bregma is slightly above the brain
    dirz = lambda-bregma;      % lambda is slightly above the brain
    diry = bregma-midsagittal; % midsagittal should be somewhere in the middle of the brain
    dirx = cross(diry,dirz);   % compute the positive x-direction, i.e. towards the right
    diry = cross(dirz,dirx);   % update the y-direction
    dirx = dirx/norm(dirx);
    diry = diry/norm(diry);
    dirz = dirz/norm(dirz);
  otherwise
    error('unrecognized headcoordinate system "%s"', coordsys);
end

% use the extra point to validate that it is a right-handed coordinate system
if ~isempty(extrapoint)
  dirq = extrapoint-origin;
  dirq = dirq/norm(dirq);
  if any(strcmp(coordsys, {'als_ctf' 'ctf' 'bti' '4d' 'yokogawa' 'ras_itab' 'itab' 'neuromag'}))
    phi = dirq(:)'*dirz(:);
    if sign(phi)<0
      warning('the input coordinate system seems left-handed, flipping z-axis to keep the transformation matrix consistent');
      dirz = -dirz;
    end
  elseif any(strcmp(coordsys, {'ras_tal' 'tal' 'spm'}))
    phi = dirq(:)'*dirx(:);
    if sign(phi)<0
      warning('the input coordinate system seems left-handed, flipping x-axis to keep the transformation matrix consistent');
      dirx = -dirx;
    end
  else
    warning('the extra input coordinate is not used with coordsys "%s"', coordsys);
  end
end

% compute the rotation matrix
rot = eye(4);
rot(1:3,1:3) = inv(eye(3) / [dirx; diry; dirz]);
% compute the translation matrix
tra = eye(4);
tra(1:4,4)   = [-origin(:); 1];
% compute the full homogeneous transformation matrix from these two
h = rot * tra;
