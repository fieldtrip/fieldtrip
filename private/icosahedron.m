function [pnt, dhk] = icosahedron();

% ICOSAHEDRON creates an icosahedron
%
% [pnt, dhk] = icosahedron
% creates an icosahedron with 12 vertices and 20 triangles
% 
% See also OCTAHEDRON, ICOSAHEDRON42, ICOSAHEDRON162, ICOSAHEDRON642, ICOSAHEDRON2562

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: icosahedron.m,v $
% Revision 1.4  2006/07/26 11:03:38  roboos
% added "see also octahedron"
%
% Revision 1.3  2003/03/11 15:35:20  roberto
% converted all files from DOS to UNIX
%
% Revision 1.2  2003/03/04 21:46:18  roberto
% added CVS log entry and synchronized all copyright labels
%

dhk = [
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

pnt = zeros(12, 3);

rho=0.4*sqrt(5);
phi=2*pi*(0:4)/5;

pnt( 1, :) = [0 0  1];			% top point

pnt(2:6, 1) = rho*cos(phi)';
pnt(2:6, 2) = rho*sin(phi)';
pnt(2:6, 3) = rho/2;

pnt(7:11, 1) = rho*cos(phi - pi/5)';
pnt(7:11, 2) = rho*sin(phi - pi/5)';
pnt(7:11, 3) = -rho/2;

pnt(12, :) = [0 0 -1];			% bottom point

