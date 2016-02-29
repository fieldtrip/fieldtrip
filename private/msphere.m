function [pnt, tri] = msphere(N)

% MSPHERE returns a triangulated sphere with approximately M vertices
% that are nicely distributed over the sphere. The vertices are aligned
% along equally spaced horizontal contours according to an algorithm of
% Dave Russel.
% 
% Use as
%  [pnt, tri] = msphere(M)
%
% See also SPHERE, NSPHERE, ICOSAHEDRON, REFINE

% Copyright (C) 1994, Dave Rusin 
% 
% The code of this function is in the public domain and originates 
% from the sci.math newsgroup. The full news message can be found
% below. This MATLAB implementation was made by Robert Oostenveld.
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

  % put a single vertex at the top
  phi = [0];
  th  = [0];

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
pnt = [x' y' z'];
tri = convhulln(pnt);

fprintf('returning %d vertices, %d triangles\n', size(pnt,1), size(tri,1));

