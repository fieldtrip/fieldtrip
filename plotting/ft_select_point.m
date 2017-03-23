function [selected] = ft_select_point(pos, varargin)

% FT_SELECT_POINT helper function for selecting a one or multiple points in the
% current figure using the mouse. It returns a list of the [x y] coordinates of the
% selected points.
%
% Use as
%   [selected] = ft_select_point(pos, ...)
%
% Optional input arguments should come in key-value pairs and can include
%   'multiple'   = true/false, make multiple selections, pressing "q" on the keyboard finalizes the selection (default = false)
%   'nearest'    = true/false (default = true)
%
% Example
%   pos = randn(10,2);
%   figure
%   plot(pos(:,1), pos(:,2), '.')
%   ft_select_point(pos)
%
% See also FT_SELECT_BOX, FT_SELECT_CHANNEL, FT_SELECT_POINT3D, FT_SELECT_RANGE, FT_SELECT_VOXEL 

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
nearest  = ft_getopt(varargin, 'nearest',   true);
multiple = ft_getopt(varargin, 'multiple',  false);

% ensure that it is boolean
nearest  = istrue(nearest);
multiple = istrue(multiple);

if multiple
  fprintf('select multiple points by clicking in the figure, press "q" if you are done\n');
end

x = [];
y = [];
done = false;

selected = zeros(0,3);

% ensure that "q" is not the current character, which happens if you reuse the same figure
set(gcf, 'CurrentCharacter', 'x')

while ~done
  k     = waitforbuttonpress;
  point = get(gca,'CurrentPoint');     % button down detected
  key   = get(gcf,'CurrentCharacter'); % which key was pressed (if any)?
  if strcmp(key, 'q')
    % we are done with the clicking
    done = true;
  else
    % add the current point
    x(end+1) = point(1,1);
    y(end+1) = point(1,2);
  end
  
  if ~multiple
    done = true;
  end
end

if nearest && ~isempty(pos)
  % determine the points that are the nearest to the displayed points
  selected = [];
  
  % compute the distance between the points to get an estimate of the tolerance
  dp = dist(pos');
  dp = triu(dp, 1);
  dp = dp(:);
  dp = dp(dp>0);
  % allow for some tolerance in the clicking
  dp = median(dp);
  tolerance = 0.3*dp;
  
  for i=1:length(x)
    % compute the distance between the clicked position and all points
    dx = pos(:,1) - x(i);
    dy = pos(:,2) - y(i);
    dd = sqrt(dx.^2 + dy.^2);
    [d, i] = min(dd);
    if d<tolerance
      selected(end+1,:) = pos(i,:);
    end
  end
  
else
  selected = [x(:) y(:)];
end % if nearest

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function serves as a replacement for the dist function in the Neural
% Networks toolbox.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d] = dist(x)
n = size(x,2);
d = zeros(n,n);
for i=1:n
  for j=(i+1):n
    d(i,j) = sqrt(sum((x(:,i)-x(:,j)).^2));
    d(j,i) = d(i,j);
  end
end

