function ft_select_range(handle, eventdata, varargin)

% FT_SELECT_RANGE is a helper function that can be used as callback function
% in a figure. It allows the user to select a horizontal or a vertical
% range, or one or multiple boxes.
%
% The callback function (and it's arguments) specified in callback is called
% on a left-click inside a selection, or using the right-click context-menu.
% The callback function will have as its first-to-last input argument the range of
% all selections. The last input argument is either empty, or, when using the context
% menu, a label of the item clicked.
% Context menus are shown as the labels presented in the input. When activated,
% the callback function is called, with the last input argument being the label of
% the selection option.
%
% Input arguments:
%   'event'       = string, event used as hook.
%   'callback'    = function handle or cell-array containing function handle and additional input arguments
%   'contextmenu' = cell-array containing labels shown in right-click menu
%   'multiple'    = boolean, allowing multiple selection boxes or not
%   'xrange'      = boolean, xrange variable or not
%   'yrange'      = boolean, yrange variable or not
%   'clear'       = boolean
%
% Example
%   x = randn(10,1);
%   y = randn(10,1);
%   figure; plot(x, y, '.');
%
% The following example allows multiple horizontal and vertical selections to be made
%   set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', true, 'callback', @disp});
%   set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', true, 'callback', @disp});
%   set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', true, 'callback', @disp});
%
% The following example allows a single horizontal selection to be made
%   set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', false, 'xrange', true, 'yrange', false, 'callback', @disp});
%   set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', false, 'xrange', true, 'yrange', false, 'callback', @disp});
%   set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', false, 'xrange', true, 'yrange', false, 'callback', @disp});
%
% The following example allows a single point to be selected
%   set(gcf, 'WindowButtonDownFcn',   {@ft_select_range, 'event', 'WindowButtonDownFcn',   'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp});
%   set(gcf, 'WindowButtonMotionFcn', {@ft_select_range, 'event', 'WindowButtonMotionFcn', 'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp});
%   set(gcf, 'WindowButtonUpFcn',     {@ft_select_range, 'event', 'WindowButtonUpFcn',     'multiple', false, 'xrange', false, 'yrange', false, 'callback', @disp});
%
% See also FT_SELECT_BOX, FT_SELECT_CHANNEL, FT_SELECT_POINT, FT_SELECT_POINT3D, FT_SELECT_VOXEL

% Copyright (C) 2009-2020, Robert Oostenveld
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

% get the optional arguments
event       = ft_getopt(varargin, 'event');
callback    = ft_getopt(varargin, 'callback');
multiple    = ft_getopt(varargin, 'multiple', false);
xrange      = ft_getopt(varargin, 'xrange',   true);
yrange      = ft_getopt(varargin, 'yrange',   true);
clear       = ft_getopt(varargin, 'clear',    false);
contextmenu = ft_getopt(varargin, 'contextmenu'); % this will be displayed following a right mouse click

% convert 'yes/no' string to boolean value
multiple  = istrue(multiple);
xrange    = istrue(xrange);
yrange    = istrue(yrange);
clear     = istrue(clear);

% get the figure handle, dependent on MATLAB version
if ft_platform_supports('graphics_objects')
  while ~isa(handle, 'matlab.ui.Figure')
    handle = get(handle, 'parent');
  end
else
  while ~isequal(handle, 0)
    handle = get(handle, 'parent');
  end
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

xLim = get(gca, 'XLim');
yLim = get(gca, 'YLim');

% limit cursor coordinates
if p(1)<xLim(1), p(1)=xLim(1); end
if p(1)>xLim(2), p(1)=xLim(2); end
if p(2)<yLim(1), p(2)=yLim(1); end
if p(2)>yLim(2), p(2)=yLim(2); end

% determine whether the user is currently making a selection
selecting = numel(userData.range)>0 && any(isnan(userData.range(end,:)));
pointonly = ~xrange && ~yrange;

if pointonly && multiple
  ft_warning('multiple selections are not possible for a point');
  multiple = false;
end

% setup contextmenu
if ~isempty(contextmenu)
  if isempty(get(handle,'uicontextmenu'))
    hcmenu    = uicontextmenu(handle);
    hcmenuopt = nan(1,numel(contextmenu));
    for icmenu = 1:numel(contextmenu)
      hcmenuopt(icmenu) = uimenu(hcmenu, 'label', contextmenu{icmenu}, 'callback', {@evalcontextcallback, callback{:}, []}); % empty matrix is placeholder, will be updated to userdata.range
    end
  end
  if ~exist('hcmenuopt','var')
    hcmenuopt = get(get(handle,'uicontextmenu'),'children'); % uimenu handles, used for switchen on/off and updating
  end
  
  % setting associations for all clickable objects
  % this out to be pretty fast, if this is still to slow in some cases, the code below has to be reused
  if ~exist('hcmenu','var')
    hcmenu = get(handle,'uicontextmenu');
  end
  set(findobj(handle,'hittest','on'), 'uicontextmenu',hcmenu);
  % to be used if above is too slow
  % associations only done once. this might be an issue in some cases, cause when a redraw is performed in the original figure (e.g. databrowser), a specific assocations are lost (lines/patches/text)
  %   set(get(handle,'children'),'uicontextmenu',hcmenu);
  %   set(findobj(handle,'type','text'), 'uicontextmenu',hcmenu);
  %   set(findobj(handle,'type','patch'),'uicontextmenu',hcmenu);
  %   set(findobj(handle,'type','line'), 'uicontextmenu',hcmenu); % fixme: add other often used object types that ft_select_range is called upon
end

% get last-used-mouse-button
lastmousebttn = get(gcf,'selectiontype');

switch lower(event)
  
  case lower('WindowButtonDownFcn')
    switch lastmousebttn
      case 'normal' % left click
        
        if inSelection(p, userData.range)
          % the user has clicked in one of the existing selections
          evalCallback(callback, userData.range);
          if clear
            delete(userData.box(ishandle(userData.box)));
            userData.range = [];
            userData.box   = [];
            set(handle, 'Pointer', 'crosshair');
            if ~isempty(contextmenu) && ~pointonly
              set(hcmenuopt,'enable', 'off')
            end
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
    end
    
    
  case lower('WindowButtonUpFcn')
    switch lastmousebttn
      case 'normal' % left click
        
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
          % update contextmenu callbacks
          if ~isempty(contextmenu)
            updateContextCallback(hcmenuopt, callback, userData.range)
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
    end
    
    
  case lower('WindowButtonMotionFcn')
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
      %set(userData.box(end), 'EraseMode', 'xor');
      set(userData.box(end), 'LineStyle', '--');
      set(userData.box(end), 'LineWidth', 1.5);
      set(userData.box(end), 'Visible', 'on');
      
    else
      % update the cursor
      if inSelection(p, userData.range)
        set(handle, 'Pointer', 'hand');
        if ~isempty(contextmenu)
          set(hcmenuopt,'enable','on')
        end
      else
        set(handle, 'Pointer', 'crosshair');
        if ~isempty(contextmenu)
          set(hcmenuopt,'enable','off')
        end
      end
    end
    
    
  otherwise
    ft_error('unexpected event "%s"', event);
    
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
% no context menu item was clicked, set to empty
cmenulab = [];

if ~isempty(callback)
  if isa(callback, 'cell')
    % the callback specifies a function and additional arguments
    funhandle = callback{1};
    funargs   = callback(2:end);
    feval(funhandle, funargs{:}, val, cmenulab);
  else
    % the callback only specifies a function
    funhandle = callback;
    feval(funhandle, val);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateContextCallback(hcmenuopt, callback, val)
if ~isempty(callback)
  if isa(callback, 'cell')
    % the callback specifies a function and additional arguments
    funhandle = callback{1};
    funargs   = callback(2:end);
    callback  = {funhandle, funargs{:}, val};
  else
    % the callback only specifies a function
    funhandle = callback;
    callback  = {funhandle, val};
  end
  for icmenu = 1:numel(hcmenuopt)
    set(hcmenuopt(icmenu),'callback',{@evalcontextcallback, callback{:}})
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evalcontextcallback(hcmenuopt, eventdata, varargin)

% delete selection box if present
% get parent (uimenu -> uicontextmenu -> parent)
parent = get(get(hcmenuopt,'parent'),'parent'); % fixme: isn't the parent handle always input provided in the callback?
userData = getappdata(parent, 'select_range_m');
if ishandle(userData.box)
  if any(~isnan([get(userData.box,'ydata') get(userData.box,'xdata')]))
    delete(userData.box(ishandle(userData.box)));
    userData.range = [];
    userData.box   = [];
    set(parent, 'Pointer', 'crosshair');
    setappdata(parent, 'select_range_m', userData);
  end
end

% get contextmenu name
cmenulab = get(hcmenuopt,'label');
if numel(varargin)>1
  % the callback specifies a function and additional arguments
  funhandle = varargin{1};
  funargs   = varargin(2:end);
  feval(funhandle, funargs{:}, cmenulab);
else
  % the callback only specifies a function
  funhandle = varargin{1};
  feval(funhandle, val, cmenulab);
end
