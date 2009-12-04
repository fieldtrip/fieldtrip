function [tra] = transfer_elec(pnt, tri, el); 

% TRANSFER_ELEC is the transfermatrix from vertex to electrode potentials
% using bi-linear interpolation over the triangles
%
% tra = transfer_elec(pnt, tri, el)
%
% the Nx3 matrix el shold contain [tri, la, mu] for each electrode
%
% See also PROJECT_ELEC

% Copyright (C) 1998-2002, Robert Oostenveld
% 
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

Npnt = size(pnt,1);
Ntri = size(tri,1);
Nel  = size(el,1);

tra = zeros(Nel, Npnt);

for i=1:Nel
  tra(i, tri(el(i,1), 1)) = 1 - el(i,2) - el(i,3);
  tra(i, tri(el(i,1), 2)) = el(i,2);
  tra(i, tri(el(i,1), 3)) = el(i,3);
end

