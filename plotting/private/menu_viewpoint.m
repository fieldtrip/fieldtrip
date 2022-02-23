function ft_viewpoint(h, coordsys)

h = getparent(h);
setappdata(h, 'coordsys', coordsys);

cm = uicontextmenu(h);
viewpoint = {'top', 'bottom', 'left', 'right', 'front', 'back'};
for i=1:length(viewpoint)
  uimenu(cm, 'Text', viewpoint{i}, 'callback', @cb_viewpoint);
end

h.ContextMenu = cm;

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