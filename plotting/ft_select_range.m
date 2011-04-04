function ft_select_range(handle, eventdata, varargin)

% FT_SELECT_RANGE is a helper function that can be used as callback function
% in a figure. It allows the user to select a horizontal or a vertical
% range, or one or multiple boxes.
%
% Example
%   x = randn(10,1);
%   y = randn(10,1);
%   figure; plot(x, y, '.');
%
%   set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'});
%   set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonUpFcn'});
%   set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonMotionFcn'});
%
%   set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp, 'event', 'WindowButtonDownFcn'});
%   set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp, 'event', 'WindowButtonUpFcn'});
%   set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp, 'event', 'WindowButtonMotionFcn'});

% Copyright (C) 2009, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% get the optional arguments
event    = keyval('event',    varargin);
callback = keyval('callback', varargin);
multiple = keyval('multiple', varargin); if isempty(multiple), multiple = false; end
xrange   = keyval('xrange',   varargin); if isempty(xrange), xrange = true; end
yrange   = keyval('yrange',   varargin); if isempty(yrange), yrange = true; end
clear    = keyval('clear',   varargin);  if isempty(clear),  clear = false; end
contextmenu = keyval('contextmenu', varargin); % this will be displayed following a right mouse click

% convert 'yes/no' string to boolean value
multiple  = istrue(multiple);
xrange    = istrue(xrange);
yrange    = istrue(yrange);

p = handle;
while ~isequal(p, 0)
  handle = p;
  p = get(handle, 'parent');
end

if ishandle(handle)
  userData = getappdata(handle, 'select_range_m');
else
  userData = [];
end

if isempty(userData)
  userData.range = []; % this is a Nx4 matrix with the selection range
  userData.box   = []; % this is a Nx1 vector with the line handle
end

p = get(gca, 'CurrentPoint');
p = p(1,1:2);

abc = axis;
xLim = abc(1:2);
yLim = abc(3:4);

% limit cursor coordinates
if p(1)<xLim(1), p(1)=xLim(1); end;
if p(1)>xLim(2), p(1)=xLim(2); end;
if p(2)<yLim(1), p(2)=yLim(1); end;
if p(2)>yLim(2), p(2)=yLim(2); end;

% determine whether the user is currently making a selection
selecting = numel(userData.range)>0 && any(isnan(userData.range(end,:)));
pointonly = ~xrange && ~yrange;

if pointonly && multiple
  warning('multiple selections are not possible for a point');
  multiple = false;
end

switch lower(event)
  case 'windowbuttondownfcn'
    if inSelection(p, userData.range)
      % the user has clicked in one of the existing selections
      evalCallback(callback, userData.range);
      if clear
        delete(userData.box(ishandle(userData.box)));
        userData.range = [];
        userData.box   = [];
        set(handle, 'Pointer', 'crosshair');
      end
      
    else
      if ~multiple
        % start with a new selection
        delete(userData.box(ishandle(userData.box)));
        userData.range = [];
        userData.box   = [];
      end
      
      % add a new selection range
      userData.range(end+1,1:4) = nan;
      userData.range(end,1) = p(1);
      userData.range(end,3) = p(2);
      
      % add a new selection box
      xData = [nan nan nan nan nan];
      yData = [nan nan nan nan nan];
      userData.box(end+1) = line(xData, yData);
    end
    
  case 'windowbuttonupfcn'
    if selecting
      % select the other corner of the box
      userData.range(end,2) = p(1);
      userData.range(end,4) = p(2);
    end
    
    if multiple && ~isempty(userData.range) && ~diff(userData.range(end,1:2)) && ~diff(userData.range(end,3:4))
      % start with a new selection
      delete(userData.box(ishandle(userData.box)));
      userData.range = [];
      userData.box   = [];
    end
    
    if ~isempty(userData.range)
      % ensure that the selection is sane
      if diff(userData.range(end,1:2))<0
        userData.range(end,1:2) = userData.range(end,[2 1]);
      end
      if diff(userData.range(end,3:4))<0
        userData.range(end,3:4) = userData.range(end,[4 3]);
      end
      if pointonly
        % only select a single point
        userData.range(end,2) = userData.range(end,1);
        userData.range(end,4) = userData.range(end,3);
      elseif ~xrange
        % only select along the y-axis
        userData.range(end,1:2) = [-inf inf];
      elseif ~yrange
        % only select along the x-axis
        userData.range(end,3:4) = [-inf inf];
      end
    end
    
    if pointonly && ~multiple
      evalCallback(callback, userData.range);
      if clear
        delete(userData.box(ishandle(userData.box)));
        userData.range = [];
        userData.box   = [];
        set(handle, 'Pointer', 'crosshair');
      end
    end
    
  case 'windowbuttonmotionfcn'
    if selecting && ~pointonly
      % update the selection box
      if xrange
        x1 = userData.range(end,1);
        x2 = p(1);
      else
        x1 = xLim(1);
        x2 = xLim(2);
      end
      if yrange
        y1 = userData.range(end,3);
        y2 = p(2);
      else
        y1 = yLim(1);
        y2 = yLim(2);
      end
      
      xData = [x1 x2 x2 x1 x1];
      yData = [y1 y1 y2 y2 y1];
      set(userData.box(end), 'xData', xData);
      set(userData.box(end), 'yData', yData);
      set(userData.box(end), 'Color', [0 0 0]);
      set(userData.box(end), 'EraseMode', 'xor');
      set(userData.box(end), 'LineStyle', '--');
      set(userData.box(end), 'LineWidth', 1.5);
      set(userData.box(end), 'Visible', 'on');
      
    else
      % update the cursor
      if inSelection(p, userData.range)
        set(handle, 'Pointer', 'hand');
      else
        set(handle, 'Pointer', 'crosshair');
      end
    end
    
  otherwise
    error('unexpected event "%s"', event);
    
end % switch event

% put the modified selections back into the figure
if ishandle(handle)
  setappdata(handle, 'select_range_m', userData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = inSelection(p, range)
if isempty(range)
  retval = false;
else
  retval = (p(1)>=range(:,1) & p(1)<=range(:,2) & p(2)>=range(:,3) & p(2)<=range(:,4));
  retval = any(retval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evalCallback(callback, val)
if ~isempty(callback)
  if isa(callback, 'cell')
    % the callback specifies a function and additional arguments
    funhandle = callback{1};
    funargs   = callback(2:end);
    feval(funhandle, val, funargs{:});
  else
    % the callback only specifies a function
    funhandle = callback;
    feval(funhandle, val);
  end
end

