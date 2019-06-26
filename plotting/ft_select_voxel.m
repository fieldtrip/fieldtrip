function ft_select_voxel(handle, eventdata, varargin)

% FT_SELECT_VOXEL is a helper function that can be used as callback function
% in a figure. It allows the user to select a voxel from a (resliced) 3-D volume.
%
% Use as
%   voxel = ft_select_voxel(h, eventdata, ...)
% The first two arguments are automatically passed by MATLAB to any
% callback function.
%
% Additional options should be specified in key-value pairs and can be
%   'callback'  = function handle to be executed after channels have been selected
%
% You can pass additional arguments to the callback function in a cell-array
% like {@function_handle,arg1,arg2}
%
% Example
%   % create a figure with a random 3-D volume
%   mri = rand(128,128,128);
%   ft_plot_slice(mri, 'location', [64 64 64], 'orientation', [1 1 1]);
%   view(120,30)
%   xlabel('x'); ylabel('y'); zlabel('z'); grid on
%   axis([0 128 0 128 0 128])
%   axis equal; axis vis3d
%   axis([0 128 0 128 0 128])
%
%   % add this function as the callback to make a single selection
%   set(gcf, 'WindowButtonDownFcn', {@ft_select_voxel, 'callback', @disp})
%
% Subsequently you can click in the figure and you'll see that the disp
% function is executed as callback and that it displays the selected
% voxel.
%
% See also FT_SELECT_BOX, FT_SELECT_CHANNEL, FT_SELECT_POINT, FT_SELECT_POINT3D, FT_SELECT_RANGE

% Copyright (C) 2010, Robert Oostenveld
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

% get optional input arguments
callback = ft_getopt(varargin, 'callback');
event    = ft_getopt(varargin, 'event');

% get the clicked position from the figure
% voxel = get(gca, 'CurrentPoint');

% the select3d function also takes care of projections
[p v vi face facei] = select3d(gca);
voxel = p';

if isempty(voxel)
  return
end

% execute the original callback with the selected channel as input argument
if ~isempty(callback)
  if isa(callback, 'cell')
    % the callback specifies a function and additional arguments
    funhandle = callback{1};
    funargs   = callback(2:end);
    feval(funhandle, voxel, funargs{:});
  else
    % the callback only specifies a function
    funhandle = callback;
    feval(funhandle, voxel);
  end
end
