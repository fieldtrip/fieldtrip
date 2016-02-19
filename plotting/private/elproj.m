function [proj] = elproj(pos, method)

% ELPROJ makes a azimuthal projection of a 3D electrode cloud
%  on a plane tangent to the sphere fitted through the electrodes
%  the projection is along the z-axis
%
%  [proj] = elproj([x, y, z], 'method');
%
% Method should be one of these:
%     'gnomic'
%     'stereographic'
%     'orthographic'
%     'inverse'
%     'polar'
%
% Imagine a plane being placed against (tangent to) a globe. If
% a light source inside the globe projects the graticule onto
% the plane the result would be a planar, or azimuthal, map
% projection. If the imaginary light is inside the globe a Gnomonic
% projection results, if the light is antipodal a Sterographic,
% and if at infinity, an Orthographic.
%
% The default projection is a polar projection (BESA like).
% An inverse projection is the opposite of the default polar projection.

% Copyright (C) 2000-2008, Robert Oostenveld
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

x = pos(:,1);
y = pos(:,2);
if size(pos, 2)==3
  z = pos(:,3);
end

if nargin<2
  method='polar';
end

if strcmp(method, 'orthographic')
  % this method compresses the lowest electrodes very much
  % electrodes on the bottom half of the sphere are folded inwards
  xp = x;
  yp = y;
  num = length(find(z<0));
  str = sprintf('%d electrodes may be folded inwards in orthographic projection\n', num);
  if num
    warning(str);
  end
  proj = [xp yp];

elseif strcmp(method, 'gnomic')
  % the lightsource is in the middle of the sphere
  % electrodes on the equator are projected at infinity
  % electrodes below the equator are not projected at all
  rad = mean(sqrt(x.^2 + y.^2 + z.^2));
  phi = cart2pol(x, y);
  th  = atan(sqrt(x.^2 + y.^2) ./ z);
  xp  = cos(phi) .* tan(th) .* rad;
  yp  = sin(phi) .* tan(th) .* rad;
  num = length(find(th==pi/2 | z<0));
  str = sprintf('removing %d electrodes from gnomic projection\n', num);
  if num
    warning(str);
  end
  xp(find(th==pi/2 | z<0)) = NaN;
  yp(find(th==pi/2 | z<0)) = NaN;
  proj = [xp yp];

elseif strcmp(method, 'stereographic')
  % the lightsource is antipodal (on the south-pole)
  rad = mean(sqrt(x.^2 + y.^2 + z.^2));
  z   = z + rad;
  phi = cart2pol(x, y);
  th  = atan(sqrt(x.^2 + y.^2) ./ z);
  xp  = cos(phi) .* tan(th) .* rad * 2;
  yp  = sin(phi) .* tan(th) .* rad * 2;
  num = length(find(th==pi/2 | z<0));
  str = sprintf('removing %d electrodes from stereographic projection\n', num);
  if num
    warning(str);
  end
  xp(find(th==pi/2 | z<0)) = NaN;
  yp(find(th==pi/2 | z<0)) = NaN;
  proj = [xp yp];

elseif strcmp(method, 'inverse')
  % compute the inverse projection of the default angular projection
  [th, r] = cart2pol(x, y);
  [xi, yi, zi] = sph2cart(th, pi/2 - r, 1);
  proj = [xi, yi, zi];

elseif strcmp(method, 'polar')
  % use default angular projection
  [az, el, r] = cart2sph(x, y, z);
  [x, y] = pol2cart(az, pi/2 - el);
  proj = [x, y];

else
  error('unsupported method (%s)', method);
end
