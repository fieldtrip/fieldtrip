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

% Copyright (C) 2006, Robert Oostenveld
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

if nargin<2
  r = rank(pos);
  switch r
    case 1
      ft_warning('points are lying on a line, cannot make triangulation');
      tri = zeros(0,3);
      return
    case 2
      method = 'delaunay';
    otherwise
      method = 'convhull';
  end
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
    if all(pos(:,3)==0)
      % this can happen with simulated electrode grids
      prj = pos(:,1:2);
    else
      % make a 2D triangulation of the projected points using delaunay
      prj = elproj(pos);
    end
    tri = delaunay(prj(:,1), prj(:,2));

  otherwise
    ft_error('unsupported method');
end



