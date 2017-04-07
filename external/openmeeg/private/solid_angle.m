function [varargout] = funname(varargin)

% SOLID_ANGLE of a planar triangle as seen from the origin
%
% The solid angle W subtended by a surface S is defined as the surface
% area W of a unit sphere covered by the surface's projection onto the
% sphere. Solid angle is measured in steradians, and the solid angle
% corresponding to all of space being subtended is 4*pi sterradians.
%
% Use:
%   [w] = solid_angle(v1, v2, v3)
% or
%   [w] = solid_angle(pnt, tri)
% where v1, v2 and v3 are the vertices of a single triangle in 3D or
% pnt and tri contain a description of a triangular mesh (this will
% compute the solid angle for each triangle)

% Copyright (C) 2003-2009, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first section contains the plain MATLAB implementation. The mex file
% is many times faster and this function is called so frequently (for
% large meshes), that the mex file should be used in all practical cases.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [w] = solid_angle(r1, r2, r3);
%
% if nargin==2
%   % reassign the input arguments
%   pnt = r1;
%   tri = r2;
%   npnt = size(pnt,1);
%   ntri = size(tri,1);
%   w    = zeros(ntri,1);
%   % compute solid angle for each triangle
%   for i=1:ntri
%     r1 = pnt(tri(i,1),:);
%     r2 = pnt(tri(i,2),:);
%     r3 = pnt(tri(i,3),:);
%     w(i) = solid_angle(r1, r2, r3);
%   end
%   return
% elseif nargin==3
%   % compute the solid angle for this triangle
%   cp23_x = r2(2) * r3(3) - r2(3) * r3(2);
%   cp23_y = r2(3) * r3(1) - r2(1) * r3(3);
%   cp23_z = r2(1) * r3(2) - r2(2) * r3(1);
%   nom = cp23_x * r1(1) + cp23_y * r1(2) + cp23_z * r1(3);
%   n1 = sqrt (r1(1) * r1(1) + r1(2) * r1(2) + r1(3) * r1(3));
%   n2 = sqrt (r2(1) * r2(1) + r2(2) * r2(2) + r2(3) * r2(3));
%   n3 = sqrt (r3(1) * r3(1) + r3(2) * r3(2) + r3(3) * r3(3));
%   ip12 = r1(1) * r2(1) + r1(2) * r2(2) + r1(3) * r2(3);
%   ip23 = r2(1) * r3(1) + r2(2) * r3(2) + r2(3) * r3(3);
%   ip13 = r1(1) * r3(1) + r1(2) * r3(2) + r1(3) * r3(3);
%   den = n1 * n2 * n3 + ip12 * n3 + ip23 * n1 + ip13 * n2;
%   if (nom == 0)
%     if (den <= 0)
%       w = nan;
%       return
%     end
%   end
%   w = 2 * atan2 (nom, den);
%   return
% else
%   error('invalid input');
% end

% compile the missing mex file on the fly
% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);

  if ispc
    mex -I. -c ../../src/geometry.c
    mex -I. -c solid_angle.c ; mex solid_angle.c solid_angle.obj geometry.obj
  else
    mex -I. -c ../../src/geometry.c
    mex -I. -c ../../src/solid_angle.c 
    mex -o solid_angle solid_angle.o geometry.o
  end

  cd(pwdir);
  success = true;

catch
  % compilation failed
  disp(lasterr);
  error('could not locate MEX file for %s', mexname);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was juist created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end

