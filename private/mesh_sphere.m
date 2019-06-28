function [pos, tri] = mesh_sphere(n, method)

% MESH_SPHERE creates spherical mesh, with approximately nvertices vertices
%
% Use as
%   [pos, tri] = mesh_sphere(numvertices, method)
%
% The input parameter 'n' specifies the (approximate) number of vertices.
% Once log4((n-2)/10) is an integer, the mesh will be based on an icosahedron.
% Once log4((n-2)/4) is an integer, the mesh will be based on a refined octahedron.
% Once log4((n-2)/2) is an integer, the mesh will be based on a refined tetrahedron.
% Otherwise, an msphere will be used. If n is empty, or undefined, a 12 vertex
% icosahedron will be returned.
%
% The input parameter 'method' defines which function to use when an refined
% icosahedron, octahedron or tetrahedron is not possible, and can be 'msphere'
% (default), or 'ksphere'.
%
% See also MESH_TETRAHEDRON, MESH_OCTAHEDRON, MESH_ICOSAHEDRON

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

if nargin<1 || isempty(n)
  n = 12;
end
assert(isscalar(n), 'number of vertices should be specified as a scalar');

r_ico   = log((n-2)./10)./log(4);
r_octa  = log((n-2)./4)./log(4);
r_tetra = log((n-2)./2)./log(4);

if nargin<2 || isempty(method)
  % default method is dependent on n
  if round(r_tetra)==r_tetra
    method = 'tetrahedron';
  elseif round(r_ico)==r_ico
    method = 'icosahedron';
  elseif round(r_octa)==r_octa
    method = 'octahedron';
  else
    method = 'msphere';
  end
end
assert(ischar(method), 'method should be specified as a string');

switch method
  case 'ksphere'
    [pos, tri] = ksphere(n);
    
  case 'msphere'
    [pos, tri] = msphere(n);
    
  case 'tetrahedron'
    [pos, tri] = mesh_tetrahedron;
    if r_tetra>0
      % perform an n-fold refinement
      for i=1:r_tetra
        [pos, tri] = refine(pos, tri, 'banks');
      end
      % scale all vertices to the unit sphere
      pos = pos ./ repmat(sqrt(sum(pos.^2,2)), 1,3);
    end
    
  case 'icosahedron'
    [pos, tri] = mesh_icosahedron;
    if r_ico>0
      % perform an n-fold refinement
      for i=1:r_ico
        [pos, tri] = refine(pos, tri, 'banks');
      end
      % scale all vertices to the unit sphere
      pos = pos ./ repmat(sqrt(sum(pos.^2,2)), 1,3);
    end
    
  case 'octahedron'
    [pos, tri] = mesh_octahedron;
    if r_octa>0
      % perform an n-fold refinement
      for i=1:r_octa
        [pos, tri] = refine(pos, tri, 'banks');
      end
      % scale all vertices to the unit sphere
      pos = pos ./ repmat(sqrt(sum(pos.^2,2)), 1,3);
    end
    
end % switch method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos, tri] = ksphere(N)

% KSPHERE returns a triangulated sphere with K vertices that are
% approximately evenly distributed over the sphere.
%
% This implements the recommended algorithm in spherical coordinates
% (theta, phi) according to "Distributing many points on a sphere"
% by E.B. Saff and A.B.J. Kuijlaars, Mathematical Intelligencer 19.1
% (1997) 5--11
%
% See also http://www.math.niu.edu/~rusin/known-math/97/spherefaq

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
pos = [x, y, z];
tri = convhulln(pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos, tri] = msphere(N)

% MSPHERE returns a triangulated sphere with approximately M vertices
% that are nicely distributed over the sphere. The vertices are aligned
% along equally spaced horizontal contours according to an algorithm of
% Dave Rusin.
%
% Copyright (C) 1994, Dave Rusin
%
% The code of this function is in the public domain and originates
% from the sci.math newsgroup. The full news message can be found
% below. This MATLAB implementation was made by Robert Oostenveld.
%
% From: rusin@washington.math.niu.edu (Dave Rusin)
% Newsgroups: sci.math
% Subject: Re: Holes on a sphere
% Date: 8 Nov 1994 17:08:12 GMT
%
% In article <39mrqm$8f1@usenet.ins.cwru.edu>,
% John Barmann <an311@cleveland.Freenet.Edu> wrote:
% >
% >I have a geomety/math/CAD probelm that i can't seem to crack.  I need
% >to space 14,000 to 15,000 holes on a sphere.  The sphere has a diameter
% >of 8.125 inches and is hollow.  The holes are to be .098 inches in
%
% "spaced equally" isn't an agreeable notion on a sphere, but I think I
% understand your intention. With so many holes to place, why not just
% arrange them along parallel circles, starting at the north pole? To place
% N holes on a sphere of radius R leaves each one surrounded by an
% area of 4 pi R^2 / N, equivalent to a disc of radius 2R/sqrt(N) around
% each. I would draw discs of this radius around each pole, then draw
% bands of width 4R/sqrt(N) around them, and space the points along
% the centers of the bands at distance of about 4R/sqrt(N) apart.
%
% (There are a lot of round-off problems to address, though it sounds
% like you're not particular as to the exact number of holes or an
% exactly uniform placement. If you do need to improve this distribution
% you could model this as a bunch of charged particles floating on the
% surface of the sphere and watch them repel each other, and take their
% resting places as the centers of your holes.)
%
% In standard spherical coordinates (phi, theta), place your first point
% at phi=0 (theta undefined there). Subsequent points are placed at
% phi = (k/M) pi   (k=1, 2, ..., M), where M is an integer close to
% (pi/4) sqrt(N) (say M=100); take points with theta coordinates
% theta= (j/Q) (2 pi) (j=1, 2, ..., Q) where  Q is an integer close to
% 2 M sin(phi). This gives two poles and M-1 bands of points having
% between 6 and  2M points on them. Taking M=100 gives 12731 points;
% perhaps M=106 or so would match your criteria better.
%
% (These angles are in radians of course.)

% due to roundoff difficulties, an exact match cannot be found
% iteratively search for the M which gives the best match
storeM    = [];
storelen  = [];
increaseM = 0;
while (1)
  
  % put a single vertex at the top% subfunction
  
  phi = 0;
  th  = 0;
  
  M = round((pi/4)*sqrt(N)) + increaseM;
  for k=1:M
    newphi = (k/M)*pi;
    Q = round(2*M*sin(newphi));
    for j=1:Q
      phi(end+1) = newphi;
      th(end+1)  = (j/Q)*2*pi;
      % in case of even number of contours
      if mod(M,2) && k>(M/2)
        th(end) = th(end) + pi/Q;
      end
    end
  end
  
  % put a single vertex at the bottom
  phi(end+1) = [pi];
  th(end+1)  = [0];
  
  % store this vertex packing
  storeM(end+1).th  = th;
  storeM(end  ).phi = phi;
  storelen(end+1) = length(phi);
  if storelen(end)>N
    break;
  else
    increaseM = increaseM+1;
    % fprintf('increasing M by %d\n', increaseM);
  end
end

% take the vertex packing that most closely matches the requirement
[dum, i] = min(abs(storelen-N));
th  = storeM(i).th;
phi = storeM(i).phi;

% convert from spherical to cartehsian coordinates
[x, y, z] = sph2cart(th, pi/2-phi, 1);
pos = [x' y' z'];
tri = convhulln(pos);
