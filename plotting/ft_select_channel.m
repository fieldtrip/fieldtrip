function ft_select_channel(handle, eventdata, varargin)

% FT_SELECT_CHANNEL is a helper function that can be used as callback function
% in a figure. It allows the user to select a channel. The channel labels
% are returned.
%
% Use as
%   label = ft_select_channel(h, eventdata, ...)
% The first two arguments are automatically passed by MATLAB to any
% callback function.
%
% Additional options should be specified in key-value pairs and can be
%   'callback'  = function handle to be executed after channels have been selected
%
% You can pass additional arguments to the callback function in a cell-array
% like {@function_handle,arg1,arg2}
%
% Example 1
%   % create a figure
%   figure
%   cfg = [];
%   cfg.channel = {'chan1', 'chan2', 'chan3', 'chan4'};
%   cfg.layout  = 'ordered';
%   lay = ft_prepare_layout(cfg);
%   ft_plot_layout(lay)
%
%   % add the required guidata
%   info       = guidata(gcf)
%   info.x     = lay.pos(:,1);
%   info.y     = lay.pos(:,2);
%   info.label = lay.label
%   guidata(gcf, info)
%
%   % add this function as the callback to make a single selection
%   set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'callback', @disp})
%
%   % or to make multiple selections
%   set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%   set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%   set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%
% Example 2 (executed from within a subplot
%   % create a figure
%   figure
%   subplot(2,2,1)
%   cfg = [];
%   cfg.channel = {'chan1', 'chan2', 'chan3', 'chan4'};
%   cfg.layout  = 'ordered';
%   lay = ft_prepare_layout(cfg);
%   ft_plot_layout(lay) 
%
%   % add the channel information to guidata under identifier linked to this axis
%   ident              = ['axh' num2str(round(sum(clock.*1e6)))]; % unique identifier for this axis
%   set(gca,'tag',ident);
%   info               = guidata(gcf);
%   info.(ident).x     = lay.pos(:, 1);
%   info.(ident).y     = lay.pos(:, 2);
%   info.(ident).label = lay.label;
%   guidata(gcf, info)
%
%   % add this function as the callback to make a single selection
%   set(gcf, 'WindowButtonDownFcn', {@ft_select_channel, 'callback', @disp})
%
%   % or to make multiple selections
%   set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%   set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%   set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%      
%
% Subsequently you can click in the figure and you'll see that the disp
% function is executed as callback and that it displays the selected
% channels.
%
% See also FT_SELECT_BOX, FT_SELECT_POINT, FT_SELECT_POINT3D, FT_SELECT_RANGE, FT_SELECT_VOXEL 

% Copyright (C) 2009, Robert Oostenveld
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
multiple = ft_getopt(varargin, 'multiple', false);
callback = ft_getopt(varargin, 'callback');

if istrue(multiple)
  % the selection is done using select_range, which will subsequently call select_channel_multiple
  set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', true, 'callback', {@select_channel_multiple, callback}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', true, 'callback', {@select_channel_multiple, callback}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', true, 'callback', {@select_channel_multiple, callback}, 'event', 'WindowButtonMotionFcn'});
else
  % the selection is done using select_channel_single
  pos = get(gca, 'CurrentPoint');
  pos = pos(1,1:2);
  select_channel_single(pos, callback)
end % if multiple


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to assist in the selection of a single channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_channel_single(pos, callback)
% If this function is called from a plotting function that supports interactive subplots,
% then the current axis, which was used for selection, was given a tag that starts with 'axh' followed by a unique ID.
% If so, then the data (e.g. channel info) that should be used is contained in a subfield of the info (guidata)
% structure, with field name being the 'tag' in question. 
% If no such tag is found, then the simple situation is assumed, where the data is contained in the 'root' field
tag = get(gca,'tag');
% in case subplots were used
if strncmp(tag,'axh',3)
  info   = guidata(gcf);
  x      = info.(tag).x(:);
  y      = info.(tag).y(:);
  label  = info.(tag).label(:);  
else % no subplots were used
  info   = guidata(gcf);
  x      = info.x(:);
  y      = info.y(:);
  label  = info.label(:);
end

% compute a tolerance measure
distance = sqrt(abs(sum([x y]'.*[x y]',1)));
distance = triu(distance, 1);
distance = distance(:);
distance = distance(distance>0);
distance = median(distance);
tolerance = 0.3*distance;

% compute the distance between the clicked point and all channels
dx = x - pos(1);
dy = y - pos(2);
dd = sqrt(dx.^2 + dy.^2);
[d, i] = min(dd);
if d<tolerance
  label = label{i};
  fprintf('channel "%s" selected\n', label);
else
  label = {};
  fprintf('no channel selected\n');
end

% execute the original callback with the selected channel as input argument
if ~isempty(callback)
  if isa(callback, 'cell')
    % the callback specifies a function and additional arguments
    funhandle = callback{1};
    funargs   = callback(2:end);
    feval(funhandle, label, funargs{:});
  else
    % the callback only specifies a function
    funhandle = callback;
    feval(funhandle, label);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to assist in the selection of multiple channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select_channel_multiple(callback,range,cmenulab) % last input is context menu label, see ft_select_range
% If this function is called from a plotting function that supports interactive subplots,
% then the current axis, which was used for selection, was given a tag that starts with 'axh' followed by a unique ID.
% If so, then the data (e.g. channel info) that should be used is contained in a subfield of the info (guidata)
% structure, with field name being the 'tag' in question. 
% If no such tag is found, then the simple situation is assumed, where the data is contained in the 'root' field
tag = get(gca,'tag');
% in case subplots were used
if strncmp(tag,'axh',3)
  info   = guidata(gcf);
  x      = info.(tag).x(:);
  y      = info.(tag).y(:);
  label  = info.(tag).label(:);  
else % no subplots were used
  info   = guidata(gcf);
  x      = info.x(:);
  y      = info.y(:);
  label  = info.label(:);
end

% determine which channels ly in the selected range
select = false(size(label));
for i=1:size(range,1)
  select = select | (x>=range(i, 1) & x<=range(i, 2) & y>=range(i, 3) & y<=range(i, 4));
end
label  = label(select);

% execute the original callback with the selected channels as input argument
if any(select)
  if ~isempty(callback)
    if isa(callback, 'cell')
      % the callback specifies a function and additional arguments
      funhandle = callback{1};
      funargs   = callback(2:end);
      feval(funhandle, label, funargs{:});
    else
      % the callback only specifies a function
      funhandle = callback;
      feval(funhandle, label);
    end
  end
end
