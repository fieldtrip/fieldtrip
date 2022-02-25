function ft_viewpoint(h, coordsys)

h = getparent(h);

if ~ishandle(h) || ~ft_platform_supports('uimenu')
  return
end

setappdata(h, 'coordsys', coordsys);

cm = uicontextmenu(h);
viewpoint = {'top', 'bottom', 'left', 'right', 'front', 'back'};
for i=1:length(viewpoint)
  uimenu(cm, 'Text', viewpoint{i}, 'callback', @cb_viewpoint);
end

try
  % this works for R2021b
  h.ContextMenu = cm;
catch
  % FIXME it would be nice to get this to work on older MATLAB versions
  % 2018b has h.UIContextMenu, but it seems to work a bit differently
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_viewpoint(h, eventdata)
h = getparent(h);
setviewpoint(h, getappdata(h, 'coordsys'), eventdata.Source.Text);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end