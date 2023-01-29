function menu_viewpoint(h, coordsys)

% MENU_VIEWPOINT adds a context menu to a 3D figure.
%
% See also MENU_FIELDTRIP

h = getparent(h);

if ~ishandle(h) || ~ft_platform_supports('uimenu')
  return
end

setappdata(h, 'coordsys', coordsys);

cm = uicontextmenu(h);
viewpoint = {'top', 'bottom', 'left', 'right', 'front', 'back'};
for i=1:length(viewpoint)
  uimenu(cm, 'Label', viewpoint{i}, 'callback', @cb_viewpoint);
end

if isa(h, 'matlab.ui.Figure') && isprop(h, 'ContextMenu')
  % this works for R2021b and up
  h.ContextMenu = cm;
else
  % FIXME it would be nice to get this to work on older MATLAB versions
  % (e.g.) 2018b has h.UIContextMenu, but it seems to work a bit
  % differently, but it's worth a shot as fallback option
  set(h, 'UIContextMenu', cm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_viewpoint(h, eventdata)
h  = getparent(h);
ax = findall(h, 'type', 'axes');
try
  str = eventdata.Source.Text;
catch
  str = eventdata.Source.Label;
end
setviewpoint(ax, getappdata(h, 'coordsys'), str);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end