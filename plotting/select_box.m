function [x, y] = select_box(handle, eventdata, varargin)

% SELECT_BOX helper function for selecting a rectangular region
% in the current figure using the mouse.
%
% Use as
%   [x, y] = select_box(...)
%
% It returns a 2-element vector x and a 2-element vector y
% with the corners of the selected region.
%
% Optional input arguments should come in key-value pairs and can include
%   'multiple'   true/false, make multiple selections by dragging, clicking
%                in one will finalize the selection (default = false)

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: select_box.m,v $
% Revision 1.6  2009/07/14 13:18:33  roboos
% updated channel selection, use select_range and two local helper functions, also support multiple selections
%
% Revision 1.5  2009/06/04 10:50:50  roboos
% changed handling of inputs
%
% Revision 1.4  2009/06/03 11:21:49  crimic
% bug fixes
%
% Revision 1.3  2009/05/29 15:56:08  roboos
% added input argument handling for 'multiple', the actual implementation does not support it yet
%
% Revision 1.2  2009/04/15 12:34:29  crimic
% added code of original select2d.m function
%
% Revision 1.1  2006/05/17 14:38:09  roboos
% new implementation

% get the optional arguments
multiple = keyval('multiple', varargin); if isempty(multiple), multiple = false; end

if multiple
  error('not yet implemented');
else
  k = waitforbuttonpress;
  point1 = get(gca,'CurrentPoint');    % button down detected
  finalRect = rbbox;                   % return figure units
  point2 = get(gca,'CurrentPoint');    % button up detected
  point1 = point1(1,1:2);              % extract x and y
  point2 = point2(1,1:2);
  x = sort([point1(1) point2(1)]);
  y = sort([point1(2) point2(2)]);
end


