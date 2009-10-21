function select_channel(handle, eventdata, varargin)

% SELECT_CHANNEL is a helper function that can be used as callback function
% in a figure. It allows the user to select a channel. The channel labels
% are returned.
%
% Use as
%   label = select_channel(h, eventdata, ...)
% The first two arguments are automatically passed by Matlab to any
% callback function.
%
% Additional options should be specified in key-value pairs and can be
%   'callback'  = function handle to be executed after channels have been selected
%
% You can pass additional arguments to the callback function in a cell-array
% like {@function_handle,arg1,arg2}
%
% Example
%   % create a figure
%   lay = prepare_layout([])
%   plot_lay(lay)
%
%   % add the required guidata
%   info       = guidata(gcf)
%   info.x     = lay.pos(:,1);
%   info.y     = lay.pos(:,2);
%   info.label = lay.label
%   guidata(gcf, info)
%
%   % add this function as the callback to make a single selection
%   set(gcf, 'WindowButtonDownFcn', {@select_channel, 'callback', @disp})
%
%   % or to make multiple selections
%   set(gcf, 'WindowButtonDownFcn',   {@select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%   set(gcf, 'WindowButtonUpFcn',     {@select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%   set(gcf, 'WindowButtonMotionFcn', {@select_channel, 'multiple', true, 'callback', @disp, 'event', 'WindowButtonDownFcn'})
%
% Subsequently you can click in the figure and you'll see that the disp
% function is executed as callback and that it displays the selected
% channels.

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: select_channel.m,v $
% Revision 1.4  2009/10/16 09:18:48  jansch
% ensure column representation in select_channel_multiple to avoid crash when
% called from multiplotER
%
% Revision 1.3  2009/07/14 13:18:33  roboos
% updated channel selection, use select_range and two local helper functions, also support multiple selections
%
% Revision 1.2  2009/05/12 18:10:43  roboos
% added handling of follow-up callback function
%
% Revision 1.1  2009/05/12 12:49:33  roboos
% new implementation of helper function that can be used as callback in a figure
%

% get optional input arguments
callback = keyval('callback', varargin);
multiple = keyval('multiple', varargin); if isempty(multiple), multiple = false; end

% convert 'yes/no' string to boolean value
multiple  = istrue(multiple);

if multiple
  % the selection is done using select_range, which will subsequently call select_channel_multiple
  set(gcf, 'WindowButtonDownFcn',   {@select_range, 'multiple', true, 'callback', {@select_channel_multiple, callback}, 'event', 'WindowButtonDownFcn'});
  set(gcf, 'WindowButtonUpFcn',     {@select_range, 'multiple', true, 'callback', {@select_channel_multiple, callback}, 'event', 'WindowButtonUpFcn'});
  set(gcf, 'WindowButtonMotionFcn', {@select_range, 'multiple', true, 'callback', {@select_channel_multiple, callback}, 'event', 'WindowButtonMotionFcn'});
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

info  = guidata(gcf);
x     = info.x;
y     = info.y;
label = info.label;

% compute a tolerance measure
distance = dist([x y]');
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
function select_channel_multiple(range, callback)

info   = guidata(gcf);
x      = info.x(:);
y      = info.y(:);
label  = info.label(:);

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
