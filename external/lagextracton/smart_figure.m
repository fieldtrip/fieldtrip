function [h] = smart_figure(name,show_now)
%   SMART_FIGURE
%
%       [H] = SMART_FIGURE(NAME,SHOW_NOW)
%
%       Returns the handler for the figure with a given name or creates
%       it if it does not exist
%
%   Copyright (c) 2009 Alexandre Gramfort. All rights reserved.
%
% % $Id: smart_figure.m 4 2009-08-15 21:10:35Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-08-15 17:10:35 -0400 (Sam, 15 aoû 2009) $
% $Revision: 4 $

me = 'SMART_FIGURE';

if nargin == 1
    show_now = true;
end

hw = get(0,'children');
hw = sort(hw);

for i = 1:length(hw),
  s = get(hw(i),'Name');
  if(strcmp(deblank(s),deblank(name))),
    h = hw(i);
    if show_now
        figure(h)
    end
    return;
  end
end

% we exited out without a match, make the window
h = figure;
set(h,'Name',name);

end %  function

