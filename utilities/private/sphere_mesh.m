function [pos, tri] = sphere_mesh(n, method)

% SPHERE_MESH creates spherical mesh, with approximately nvertices vertices
%
% Use as
%   [pos, tri] = sphere_mesh(numvertices, method)
%
% Where the input parameter  n specifies the (approximate) number of vertices.
% Once log4((n-2)/10) is an integer, the mesh will be based on a refined
% icosahedron, using Robert's old icosahedronXXX functionality. Otherwise,
% an msphere will be used. If n is empty, or undefined, a 12 vertex
% icosahedron will be returned. The method parameter defines which function
% to use when an refined icosahedron is not possible, and can be 'msphere'
% (default), or 'ksphere'.
%
% See also TETRAHEDRON, OCTAHEDRON

% Copyright (C) 2002, Robert Oostenveld
% Copyright (C) 2019, Robert Oostenveld and Jan-Mathijs Schoffelen
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

if nargin==0
  n = 12;
end

if nargin==1
  method = 'msphere';
end

r = log((n-2)./10)./log(4);

if round(r)==r
  
  tri = [
    1   2   3
    1   3   4
    1   4   5
    1   5   6
    1   6   2
    2   8   3
    3   9   4
    4  10   5
    5  11   6
    6   7   2
    7   8   2
    8   9   3
    9  10   4
    10  11   5
    11   7   6
    12   8   7
    12   9   8
    12  10   9
    12  11  10
    12   7  11
    ];
  
  pos = zeros(12, 3);
  
  rho=0.4*sqrt(5);
  phi=2*pi*(0:4)/5;
  
  pos( 1, :) = [0 0  1];          % top point
  
  pos(2:6, 1) = rho*cos(phi)';
  pos(2:6, 2) = rho*sin(phi)';
  pos(2:6, 3) = rho/2;
  
  pos(7:11, 1) = rho*cos(phi - pi/5)';
  pos(7:11, 2) = rho*sin(phi - pi/5)';
  pos(7:11, 3) = -rho/2;
  
  pos(12, :) = [0 0 -1];          % bottom point
  
  if r>0
    % perform an n-fold refinement
    for i=1:r
      [pos, tri] = refine(pos, tri);
      % scale all vertices to the unit sphere
      pos = pos ./ repmat(sqrt(sum(pos.^2,2)), 1, 3);
    end
  end
else
  switch method
    case 'msphere'
      [pos, tri] = msphere(n);
    case 'ksphere'
      [pos, tri] = ksphere(n);
  end
end
