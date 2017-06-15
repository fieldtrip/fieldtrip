function select3dtool(arg)
%SELECT3DTOOL A simple tool for interactively obtaining 3-D coordinates 
%
% SELECT3DTOOL(FIG) Specify figure handle
%
% Example:
%   surf(peaks);
%   select3dtool;
%   % click on surface

if nargin<1
   arg = gcf;
end

if ~ishandle(arg)
   feval(arg);
   return;
end

%% initialize gui %%
fig = arg;
figure(fig);

uistate = uiclearmode(fig);
[tool, htext] = createUI;
hmarker1 = line('marker','o','markersize',10,'markerfacecolor','k','erasemode','xor','visible','off');
hmarker2 = line('marker','o','markersize',10,'markerfacecolor','r','erasemode','xor','visible','off');

state.uistate = uistate;
state.text = htext;
state.tool = tool;
state.fig = fig;
state.marker1 = hmarker1;
state.marker2 = hmarker2;
setappdata(fig,'select3dtool',state);
setappdata(state.tool,'select3dhost',fig);

set(fig,'windowbuttondownfcn','select3dtool(''click'')');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function off

state = getappdata(gcbf,'select3dtool');

if ~isempty(state)
    delete(state.tool);
end

fig = getappdata(gcbf,'select3dhost');

if ~isempty(fig) && ishandle(fig)
    state = getappdata(fig,'select3dtool');     
    uirestore(state.uistate);
    delete(state.marker1);
    delete(state.marker2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function click

[p, v, vi] = select3d;
state = getappdata(gcbf,'select3dtool');

if ~ishandle(state.text)
    state.text = createUI;
end

if ~ishandle(state.marker1)
    state.marker1 = [];
end

if ~ishandle(state.marker2)
    state.marker2 = [];
end

setappdata(state.fig,'select3dtool',state);

if isempty(v)
    v = [nan nan nan];
    vi = nan;
    set(state.marker2,'visible','off');
else
    set(state.marker2,'visible','on','xdata',v(1),'ydata',v(2),'zdata',v(3));
end

if isempty(p)
    p = [nan nan nan];
    set(state.marker1,'visible','off');
else
    set(state.marker1,'visible','on','xdata',p(1),'ydata',p(2),'zdata',p(3));
end

% Update tool and markers
set(state.text,'string',createString(p(1),p(2),p(3),v(1),v(2),v(3),vi));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig, h] = createUI

pos = [200 200 200 200];

% Create selection tool %
fig = figure('handlevisibility','off','menubar','none','resize','off',...
    'numbertitle','off','name','Select 3-D Tool','position',pos,'deletefcn','select3dtool(''off'')');

h = uicontrol('style','text','parent',fig,'string',createString(0,0,0,0,0,0,0),...
    'units','norm','position',[0 0 1 1],'horizontalalignment','left');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [str] = createString(px,py,pz,vx,vy,vz,vi)

str = sprintf('  Position:\n  X  %f\n  Y:  %f\n  Z:  %f  \n\n  Vertex:\n  X:  %f\n  Y:  %f\n  Z:  %f  \n\n  Vertex Index:\n  %d',px,py,pz,vx,vy,vz,vi);
