function [tri] = projecttri(pnt, method)

% PROJECTTRI makes a closed triangulation of a list of vertices by
% projecting them onto a unit sphere and subsequently by constructing
% a convex hull triangulation.
%
% Use as
%   [tri] = projecttri(pnt, method)
% The optional method argument can be 'convhull' (default) or 'delaunay'.

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
  method = 'convhull';
end

switch method
  case 'convhull'
    ori = (min(pnt) + max(pnt))./2;
    pnt(:,1) = pnt(:,1) - ori(1);
    pnt(:,2) = pnt(:,2) - ori(2);
    pnt(:,3) = pnt(:,3) - ori(3);
    nrm = sqrt(sum(pnt.^2, 2));
    pnt(:,1) = pnt(:,1)./nrm;
    pnt(:,2) = pnt(:,2)./nrm;
    pnt(:,3) = pnt(:,3)./nrm;
    tri = convhulln(pnt);
  case 'delaunay'
    % make a 2d triangulation of the projected points using delaunay
    prj = elproj(pnt);
    tri = delaunay(prj(:,1), prj(:,2));
  otherwise
    error('unsupported method');
end



