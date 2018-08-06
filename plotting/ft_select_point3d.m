function [selected] = ft_select_point3d(bnd, varargin)

% FT_SELECT_POINT3D helper function for selecting one or multiple points on a 3D mesh
% using the mouse. It returns a list of the [x y z] coordinates of the selected
% points.
%
% Use as
%   [selected] = ft_select_point3d(bnd, ...)
%
% Optional input arguments should come in key-value pairs and can include
%   'multiple'    = true/false, make multiple selections, pressing "q" on the keyboard finalizes the selection (default = false)
%   'nearest'     = true/false (default = true)
%   'marker'      = character or empty, for example '.', 'o' or 'x' (default = [])
%   'markersize'  = scalar, the size of the marker (default = 10)
%   'markercolor' = character, for example 'r', 'b' or 'g' (default = 'k')
%
% Example
%   [pos, tri] = icosahedron162;
%   bnd.pos = pos;
%   bnd.tri = tri;
%   ft_plot_mesh(bnd)
%   camlight
%   ... do something here
%
% See also FT_SELECT_BOX, FT_SELECT_CHANNEL, FT_SELECT_POINT, FT_SELECT_RANGE, FT_SELECT_VOXEL

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

% get optional input arguments
nearest     = ft_getopt(varargin, 'nearest',  true);
multiple    = ft_getopt(varargin, 'multiple', false);
marker      = ft_getopt(varargin, 'marker', []);
markersize  = ft_getopt(varargin, 'markersize', 10);
markercolor  = ft_getopt(varargin, 'markercolor', 'k');

% ensure that it is boolean
nearest  = istrue(nearest);
multiple = istrue(multiple);

% get the object handles
h = get(gca, 'children');

% select the correct objects
iscorrect = false(size(h));
for i=1:length(h)
  try
    pos = get(h(i),'vertices');
    tri = get(h(i),'faces');
    if ~isempty(bnd) && isequal(bnd.pos, pos) && isequal(bnd.tri, tri)
      % it is the same object that the user has plotted before
      iscorrect(i) = true;
    elseif isempty(bnd)
      % assume that it is the same object that the user has plotted before
      iscorrect(i) = true;
    end
  end
end
h = h(iscorrect);

if isempty(h) && ~isempty(bnd)
  figure
  ft_plot_mesh(bnd);
  camlight
  selected = ft_select_point3d(bnd, varargin{:});
  return
end

if length(h)>1
  ft_warning('using the first patch object in the figure');
  h = h(1);
end

selected = zeros(0,3);


% everything is added to the current figure
holdflag = ishold;
if ~holdflag
  hold on
end

done = false;
az = 0;
el = 0;
view(az,el);

while ~done
  k = waitforbuttonpress;
  [p, v, vi, facev, facei] = select3d(h);
  if k == 1 %checks if waitforbuttonpress was a key
    key = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?
    if strcmp(key, 'q')
      % finished selecting points
      done = true;
    elseif strcmp(key, 'r')
      % remove last point
      if ~isempty(selected)
        if ~isempty(marker)
          delete(findobj('marker', '*'));
          hs = plot3(selected(1:end-1,1), selected(1:end-1,2), selected(1:end-1,3), [markercolor marker]);
          set(hs, 'MarkerSize', markersize);
        end
        fprintf('removed latest point at [%f %f %f]\n', selected(end,1), selected(end,2), selected(end,3));
        selected = selected(1:end-1,:);
      end
    elseif strcmp(key,'+')
      zoom(1.1)
    elseif strcmp(key,'-')
      zoom(0.9)
    elseif strcmp(key,'w')
      az = az+6;
      view(az,el)
    elseif strcmp(key,'a')
      el = el+6;
      view(az,el)
    elseif strcmp(key,'s')
      az = az-6;
      view(az,el)
    elseif strcmp(key,'d')
      el = el-6;
      view(az,el)
    end
    
  elseif ~isempty(p)
    % a new point was selected
    if nearest
      selected(end+1,:) = v;
    else
      selected(end+1,:) = p;
    end % if nearest
    
    if multiple
      fprintf('selected points at\n');
      disp(selected);
    else
      fprintf('selected point at [%f %f %f]\n', selected(end,1), selected(end,2), selected(end,3));
    end
    
    if ~isempty(marker)
      hs = plot3(selected(end,1), selected(end,2), selected(end,3), [markercolor marker]);
      set(hs, 'MarkerSize', markersize);
    end
  end
  
  if ~multiple
    done = true;
  end
end

if ~holdflag
  hold off
end

