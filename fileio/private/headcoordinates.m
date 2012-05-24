function [h, flag] = headcoordinates(nas, lpa, rpa, extrapoint, flag)

% HEADCOORDINATES returns the homogeneous coordinate transformation matrix
% that converts the specified fiducials in any coordinate system (e.g. MRI)
% into the rotated and translated headcoordinate system.
%
% [h, coordsys] = headcoordinates(nas, lpa, rpa,    flag) or
% [h, coordsys] = headcoordinates(pt1, pt2, pt3,    flag)
% [h, coordsys] = headcoordinates(ac,  pc,  zpoint, flag)
%
% [h, coordsys] = headcoordinates(nas, lpa, rpa, extrapoint, flag)
% 
% [h, coordsys] = headcoordinates(nas, lpa, rpa, isrighthanded, flag)
%
% The optional flag determines how the direction of the axes and the
% location of the origin will be specified
%   according to CTF conventions:             flag = 'ctf' (default)
%   according to 4D-neuroimaging conventions: flag = '4d' or 'bti'
%   according to YOKOGAWA conventions:        flag = 'yokogawa'
%   according to ASA conventions:             flag = 'asa'
%   according to NEUROMAG conventions:        flag = 'itab'
%   according to ITAB conventions:            flag = 'neuromag'
%   according to FTG conventions:             flag = 'ftg'
%   according to Talairach conventions:       flag = 'tal'
%   according to SPM conventions:             flag = 'spm'
%
% If the function is provided with 5 input arguments, the fourth argument
% being a 3 element vector, this extrapoint will be assumed to have a
% positive Z-coordinate, and will be used to ensure correct orientation of
% the z-axis (ctf, 4d, yokogawa, itab, neuromag) or X-axis (tal, spm) NOT
% YET IMPLEMENTED. This could result in the handedness of the
% transformation to be changed, but ensures consistency with the handedness
% of the input coordinate system.
%
% If the function is provided with 5 input arguments, the fourth argument
% being a boolean, determines whether the output coordinate axes are
% right-handed (true) or not (false).
%
% The CTF/4D/YOKOGAWA coordinate system defined as follows:
%   the origin is exactly between lpa and rpa
%   the X-axis goes towards nas
%   the Y-axis goes approximately towards lpa, orthogonal to X and in the plane spanned by the fiducials
%   the Z-axis goes approximately towards the vertex, orthogonal to X and Y
%
% The ASA coordinate system is defined as follows:
%   the origin is at the orthogonal intersection of the line from rpa-lpa and the line trough nas
%   the X-axis goes towards nas
%   the Y-axis goes through rpa and lpa
%   the Z-axis goes approximately towards the vertex, orthogonal to X and Y
%
% The ITAB/NEUROMAG coordinate system is defined as follows:
%   the X-axis is from the origin towards the RPA point (exactly through)
%   the Y-axis is from the origin towards the nasion (exactly through)
%   the Z-axis is from the origin upwards orthogonal to the XY-plane
%   the origin is the intersection of the line through LPA and RPA and a line orthogonal to L passing through the nasion
%
% The FTG headcoordinate system is defined as:
%   the origin corresponds with pt1
%   the x-axis is along the line from pt1 to pt2
%   the z-axis is orthogonal to the plane spanned by pt1, pt2 and pt3
%
% The Talairach/SPM headcoordinate system is defined as:
%   the origin corresponds with the anterior commissure
%   the Y-axis is along the line from the posterior commissure to the anterior commissure
%   the Z-axis is towards the vertex, in between the hemispheres
%   the X-axis is orthogonal to the midsagittal-plane, positive to the right
%
% See also FT_ELECTRODEREALIGN, FT_VOLUMEREALIGN

% Copyright (C) 2003-2011 Robert Oostenveld
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

% check whether function call is old or new style
if nargin<4
  flag       = 0;
  extrapoint = true;
elseif nargin==4
  % old style
  flag       = extrapoint;
  extrapoint = true;
elseif nargin==5
  % do nothing
end

if isnumeric(flag)
  % these are for backward compatibility, but should preferably not be used any more
  if flag==0,
    flag = 'ctf';
  elseif flag==1,
    flag = 'asa';
  elseif flag==2,
    flag = 'ftg';
  else
    error('if flag is numeric, it should assume one of the values 0/1/2');
  end
end

% ensure that they are row vectors
lpa = lpa(:)';
rpa = rpa(:)';
nas = nas(:)';

% compute the origin and direction of the coordinate axes in MRI coordinates
switch lower(flag)
  case {'als_ctf' 'ctf' 'bti' '4d' 'yokogawa'}
    % follow CTF convention
    origin = (lpa+rpa)/2;
    dirx = nas-origin;
    dirx = dirx/norm(dirx);
    dirz = cross(dirx,lpa-rpa);
    dirz = dirz/norm(dirz);
    diry = cross(dirz,dirx);
  case 'als_asa'
    % follow ASA convention
    dirz = cross(nas-rpa, lpa-rpa);
    diry = lpa-rpa;
    dirx = cross(diry,dirz);
    dirz = dirz/norm(dirz);
    diry = diry/norm(diry);
    dirx = dirx/norm(dirx);
    origin = rpa + dot(nas-rpa,diry)*diry;
  case {'ras_itab' 'itab' 'neuromag'}
    dirz = cross(rpa-lpa,nas-lpa);
    dirx = rpa-lpa;
    diry = cross(dirz,dirx);
    dirz = dirz/norm(dirz);
    diry = diry/norm(diry);
    dirx = dirx/norm(dirx);
    origin = lpa + dot(nas-lpa,dirx)*dirx;
  case 'ftg'
    % rename the marker points for convenience
    pt1 = nas; pt2 = lpa; pt3 = rpa;
    clear nas lpa rpa
    % follow FTG conventions
    origin = pt1;
    dirx = pt2-origin;
    dirx = dirx/norm(dirx);
    diry = pt3-origin;
    dirz = cross(dirx,diry);
    dirz = dirz/norm(dirz);
    diry = cross(dirz,dirx);
  case {'ras_tal' 'tal' 'spm'}
    % rename the marker points for convenience
    ac = nas; pc = lpa; xzpoint = rpa;
    clear nas lpa rpa
    origin = ac;
    diry   = ac-pc;
    diry   = diry/norm(diry);
    dirz   = xzpoint-ac;
    dirx   = cross(diry,dirz);
    dirx   = dirx/norm(dirx);
    dirz   = cross(dirx,diry);
  otherwise
    error('unrecognized headcoordinate system requested');
end

if numel(extrapoint)==3
  dirq = extrapoint-origin;
  dirq = dirq/norm(dirq);
  if any(strcmp(lower(flag), {'als_ctf' 'ctf' 'bti' '4d' 'yokogawa' 'ras_itab' 'itab' 'neuromag'}))
    phi = dirq(:)'*dirz(:);
    if sign(phi)<0
      warning('the input coordinate system seems left-handed, flipping z-axis to keep the transformation matrix consistent');
      dirz = -dirz;
    end
  elseif any(strcmp(lower(flag), {'ras_tal' 'tal' 'spm'}))
    phi = dirq(:)'*dirx(:);
    if sign(phi)<0
      warning('the input coordinate system seems left-handed, flipping x-axis to keep the transformation matrix consistent');
    end
  else
    warning('the extra input coordinate is not used');
  end
elseif ~istrue(extrapoint)
  % left-handed axes system is requested
  warning('flipping the direction of the z-axis to create left-handed axes');
  dirz = -dirz;
end

% compute the rotation matrix
rot = eye(4);
rot(1:3,1:3) = inv(eye(3) / [dirx; diry; dirz]);
% compute the translation matrix
tra = eye(4);
tra(1:4,4)   = [-origin(:); 1];
% compute the full homogeneous transformation matrix from these two
h = rot * tra;
