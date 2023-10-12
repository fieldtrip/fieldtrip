function [indx, mind] = knnsearch(pos1, pos2)

% KNNSEARCH finds the nearest neighbor in X for each point in Y.
%
% This function serves as a drop-in replacement for the knnsearch function in the
% stats toolbox.
%
% Use as
%   indx = knnsearch(p1, p2)
% where p1 and p2 are Nx3 and Mx3 matrices with vertices in 3D space. This returns a
% Mx1 vector with indices that point into x.
%
% the idea is that the distance between 2 points is:
%
% sqrt(sum((p1(x,y,z)-p2(x,y,z)).^2)
%
% since we are dealing with relative distances, we can get rid of the sqrt:
% so we need to compute:
%
% sum((p1(x,y,z)-p2(x,y,z)).^2)
%
% this is the same as:
%
% (p1x-p2x)^2 + (p1y-p2y)^2 + (p1z-p2z)^2
%
% or, equivalently:
%
% p1x^2 + p2x^2 - 2*p1x*p2x+ ...
%
% reordering:
%
% (p1x^2 + p1y^2 + p1z^2) + cross-terms + (p2x^2 + p2y^2 + p2z^2)
%
% the last term between brackets is the same for each position-of-interest:
% so it does not change the relative distance, and the first term between
% brackets only needs to be computed once (below denoted as the 'offset'
% variable.
%
% See also DIST, DSEARCHN

if nargin~=2
  error('Incorrect number of input arguments.');
end

n1 = size(pos1,1);
n2 = size(pos2,1);

% transpose once, to speed up matrix computations
pos2 = pos2';

indx   = zeros(n2,1);
mind   = inf(1,n2);
offset = sum(pos1.^2,2);

% compute up to 1e6 pairwise distances at any given time
% this is needed to keep the memory within bounds
chunksize = round(1e6/n2);
chunks    = [(0:chunksize:(n1-1)) n1];

% loop across blocks of headmodel points, and iteratively update the
% index to the nearest sourcemodel point, based on the shortcut heuristic
% explained above

for k = 1:(numel(chunks)-1)
  iy = (chunks(k)+1):chunks(k+1);

  thisd  = offset(iy)./2 - pos1(iy,:)*pos2;
  [m, i] = min(thisd, [], 1);
  issmaller = m<mind;
  mind(issmaller) = m(issmaller);
  indx(issmaller) = iy(i(issmaller));
end

mind = mind';
