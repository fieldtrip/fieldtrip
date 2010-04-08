function uidisplaytext(str, title)

% UIDISPLAYTEXT opens a figure for displaying multi-line text
% in an "edit" user interface control element.
%
% Use as
%   uidisplaytext(str, title)

% Copyright (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if nargin<2
  title = 'unknown';
end

f = figure;
h = uicontrol('style', 'edit');

set(f, 'toolBar', 'none')
set(f, 'menuBar', 'none')

set(f, 'Name', title)
set(f, 'NumberTitle', 'off')

set(h, 'Units', 'normalized');
set(h, 'Position', [0 0 1 1]);
set(h, 'backgroundColor', [1 1 1]);
set(h, 'HorizontalAlign', 'left');
set(h, 'max', 2);
set(h, 'min', 0);
set(h, 'FontName', 'Courier');
set(h, 'FontSize', 12);
set(h, 'string', str);
