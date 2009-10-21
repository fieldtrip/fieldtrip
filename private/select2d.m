function [x, y] = select2d;

% SELECT2D helper function for selecting a rectangular region
% in the current figure using the mouse.
%
% Use as
%   [x, y] = select2d
%
% It returns a 2-element vector x and a 2-element vector y 
% with the corners of the selected region.

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: select2d.m,v $
% Revision 1.1  2006/05/17 14:38:09  roboos
% new implementation
%

k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
x = sort([point1(1) point2(1)]);
y = sort([point1(2) point2(2)]);

