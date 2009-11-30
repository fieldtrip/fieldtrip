function [h] = crosshair(pos, varargin)

% CROSSHAIR adds a crosshair at position (x,y) to the current plot
% additional options are passed to the builtin line function
% the handles of the lines are returned
% 
% h = crosshair([x,y])
% 
% see also: LINE, TEXT 

% Copyright (C) 2003, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

border = axis;

x = [ border(1) pos(1)
      border(2) pos(1) ];

y = [ pos(2)    border(3)
      pos(2)    border(4) ];

if (~ishold)
 hold on
 h = line(x, y, varargin{:});
 hold off
else
 h = line(x, y, varargin{:});
end

% make both lines the same color
set(h(2), 'Color', get(h(1), 'Color'));

