function [tri] = projecttri(pos, method)

% PROJECTTRI makes a closed triangulation of a list of vertices by
% projecting them onto a unit sphere and subsequently by constructing
% a convex hull triangulation.
%
% Use as
%   [tri] = projecttri(pos, method)
% where method is either 'convhull' (default) or 'delaunay'.
%
% See also NORMALS, PCNORMALS

% Copyright (C) 2006-2019, Robert Oostenveld
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

tmp = pos;
tmp(:,1) = tmp(:,1) - mean(tmp(:,1));
tmp(:,2) = tmp(:,2) - mean(tmp(:,2));
tmp(:,3) = tmp(:,3) - mean(tmp(:,3));
r = rank(tmp);
switch r
  case 0
    ft_warning('vertices are lying on a single point, cannot make triangulation');
    tri = zeros(0,3);
    return
  case 1
    ft_warning('vertices are lying on a straight line, cannot make triangulation');
    tri = zeros(0,3);
    return
  case 2
    if nargin<2
      method = 'delaunay';
    end
  case 3
    if nargin<2
      method = 'convhull';
    end
  otherwise
    ft_error('unexpected input');
end

switch method
  case 'convhull'
    ori = (min(pos) + max(pos))./2;
    pos(:,1) = pos(:,1) - ori(1);
    pos(:,2) = pos(:,2) - ori(2);
    pos(:,3) = pos(:,3) - ori(3);
    nrm = sqrt(sum(pos.^2, 2));
    pos(:,1) = pos(:,1)./nrm;
    pos(:,2) = pos(:,2)./nrm;
    pos(:,3) = pos(:,3)./nrm;
    tri = convhulln(pos);
    if surfaceorientation(pos, tri)<0
      % make the surface outward oriented
      tri = fliplr(tri);
    end
    
  case 'delaunay'
    if all(pos(:,1)==0)
      % this can happen with simulated electrode grids
      prj = pos(:,[2 3]);
    elseif all(pos(:,2)==0)
      % this can happen with simulated electrode grids
      prj = pos(:,[1 3]);
    elseif all(pos(:,3)==0)
      % this can happen with simulated electrode grids
      prj = pos(:,[1 2]);
    else
      % make a 2D triangulation of the projected points using delaunay
      prj = elproj(pos);
    end
    tri = delaunay(prj(:,1), prj(:,2));
    
  otherwise
    ft_error('unsupported method');
end
