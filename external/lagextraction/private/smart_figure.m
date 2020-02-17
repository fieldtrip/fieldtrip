function h = smart_figure(name)
%
%   Returns index of the figure with a given name.
%   It returns a new figure when no figure with that name exists
%   

% $Id: smart_figure.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

hw = get(0,'children');
hw = sort(hw);

for i = 1:length(hw),
  s = get(hw(i),'Name');
  if(strcmp(deblank(s),deblank(name))),
    h = hw(i);
    figure(h)
    return;
  end
end

% we exited out without a match, make the window
h = figure;
set(h,'Name',name);

return

