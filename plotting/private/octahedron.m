function [pos, tri] = octahedron(varargin)

% OCTAHEDRON
%
% Use as
%   [pos tri] = octahedron;
%
% See also TETRAHEDRON ICOSAHEDRON


pos = [
  0  0 +1
  +1  0  0
  0 +1  0
  -1  0  0
  0 -1  0
  0  0 -1
  ];

tri = [
  1 2 3
  1 3 4
  1 4 5
  1 5 2
  6 3 2
  6 4 3
  6 5 4
  6 2 5
  ];

if nargin>0
  n = varargin{1};
  % perform an n-fold refinement
  for i=1:n
    [pos, tri] = refine(pos, tri, 'banks');
  end
  % scale all vertices to the unit sphere
  pos = pos ./ repmat(sqrt(sum(pos.^2,2)), 1,3);
end