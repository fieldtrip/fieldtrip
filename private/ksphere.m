function [pnt, tri] = ksphere(N)

% KSPHERE returns a triangulated sphere with K vertices that are
% approximately evenly distributed over the sphere. 
%
% Use as
%   [pnt, tri] = ksphere(K);
%
% This implements the recommended algorithm in spherical coordinates
% (theta, phi) according to "Distributing many points on a sphere"
% by E.B. Saff and A.B.J. Kuijlaars, Mathematical Intelligencer 19.1
% (1997) 5--11
%
% See also http://www.math.niu.edu/~rusin/known-math/97/spherefaq

% Copyright (C) 2005, Robert Oostenveld
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

for k=1:N
  h = -1 + 2*(k-1)/(N-1);
  theta(k) = acos(h);
  if k==1 || k==N
    phi(k) = 0;
  else
    phi(k) = mod((phi(k-1) + 3.6/sqrt(N*(1-h^2))), (2*pi));
  end
end
az = phi;
el = theta - pi/2;
[x, y, z] = sph2cart(az', el', 1);
pnt = [x, y, z];
tri = convhulln(pnt);

