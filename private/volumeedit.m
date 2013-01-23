function [dataout] = volumeedit(data, varargin)

% VOLUMEEDIT allows for editing of a (booleanized) volume, in order to
% remove unwanted voxels. Interaction proceeds with the keyboard and the
% mouse. 

% Copyright (C) 2013, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id:$

revision = '$Id: ft_sourceplot.m 7192 2012-12-13 22:32:56Z roboos $';

datain = data;
data   = data~=0;

dim  = size(data);

xi = round(dim(1)/2);
yi = round(dim(2)/2);
zi = round(dim(3)/2);

xi = max(xi, 1); xi = min(xi, dim(1));
yi = max(yi, 1); yi = min(yi, dim(2));
zi = max(zi, 1); zi = min(zi, dim(3));

% enforce the size of the subplots to be isotropic
xdim = dim(1) + dim(2);
ydim = dim(2) + dim(3);

xsize(1) = 0.82*dim(1)/xdim;
xsize(2) = 0.82*dim(2)/xdim;
ysize(1) = 0.82*dim(3)/ydim;
ysize(2) = 0.82*dim(2)/ydim;

% create figure
h = figure;
set(h, 'color', [1 1 1]);
set(h, 'visible', 'on');
set(h, 'windowbuttondownfcn', @cb_buttonpress); 
set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
set(h, 'windowkeypressfcn',   @cb_keyboard);

% axis handles
h1 = axes('position',[0.07 0.07+ysize(2)+0.05 xsize(1) ysize(1)]);
h2 = axes('position',[0.07+xsize(1)+0.05 0.07+ysize(2)+0.05 xsize(2) ysize(1)]);
h3 = axes('position',[0.07 0.07 xsize(1) ysize(2)]);

% slice handles
hs1 = imagesc(squeeze(data(xi,:,:)),'parent',h1); colormap gray;
hs2 = imagesc(data(:,:,zi)',        'parent',h2); colormap gray;
hs3 = imagesc(squeeze(data(:,yi,:)),'parent',h3); colormap gray;
set(h1, 'tag', 'jk', 'clim', [0 1]);
set(h2, 'tag', 'ji', 'clim', [0 1]);
set(h3, 'tag', 'ik', 'clim', [0 1]);

% crosshair handles
hch1 = crosshair([zi yi], 'parent', h1, 'color', 'y');
hch2 = crosshair([xi yi], 'parent', h2, 'color', 'y');
hch3 = crosshair([zi xi], 'parent', h3, 'color', 'y');

% erasercontour
he1(1,:) = line(zi-3.5+[0 0 7 7 0],yi-3.5+[0 7 7 0 0],'color','r','parent',h1); 
he2(1,:) = line(xi-3.5+[0 0 7 7 0],yi-3.5+[0 7 7 0 0],'color','r','parent',h2);
he3(1,:) = line(zi-3.5+[0 0 7 7 0],xi-3.5+[0 7 7 0 0],'color','r','parent',h3);

% create structure to be passed to gui
opt.data          = data~=0;
opt.handlesaxes   = [h1 h2 h3];
opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
opt.handlesslice  = [hs1 hs2 hs3];
opt.handleseraser = [he1(:)';he2(:)';he3(:)'];
opt.ijk           = [xi yi zi];
opt.dim           = dim;
opt.quit          = 0;
opt.mask          = true(dim);
opt.radius        = 3;

setappdata(h, 'opt', opt);
cb_redraw(h);
    
while opt.quit==0
  uiwait(h);
  opt = getappdata(h, 'opt'); % needed to update the opt.quit
end
opt = getappdata(h, 'opt');
delete(h);

dataout = datain;
dataout(opt.mask==0) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parseKeyboardEvent(eventdata);
end
% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h,       'currentaxes');
tag     = get(curr_ax, 'tag');
switch tag
  case 'jk'
    xy = [2 3];
  case 'ji'
    xy = [2 1];
  case 'ik'
    xy = [1 3];
  otherwise
end

switch key
  case 'leftarrow'
    opt.ijk(xy(2)) = opt.ijk(xy(2)) - 1;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
  case 'rightarrow'
    opt.ijk(xy(2)) = opt.ijk(xy(2)) + 1;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
  case 'uparrow'
    opt.ijk(xy(1)) = opt.ijk(xy(1)) - 1;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
  case 'downarrow'
    opt.ijk(xy(1)) = opt.ijk(xy(1)) + 1;
    setappdata(h, 'opt', opt);
    cb_redraw(h);
  case 'd'
    % delete current voxel
    cb_eraser(h);
    cb_redraw(h);
  case 'q'
    setappdata(h, 'opt', opt);
    cb_cleanup(h);
  case 'r'
    % select the radius of the eraser box
    response = inputdlg(sprintf('radius of eraser box (in voxels)'), 'specify', 1, {num2str(opt.radius)});
    if ~isempty(response)
      response   = str2double(tokenize(response{1},' '));
      opt.radius = round(response);
      opt.radius = min(opt.radius, 100);
      opt.radius = max(opt.radius, 1);
      setappdata(h, 'opt', opt);
      cb_erasercontour(h);
    end
  case 'control+control'
    % do nothing
  case 'shift+shift'
    % do nothing
  case 'alt+alt'
    % do nothing
  otherwise
    setappdata(h, 'opt', opt);
    %cb_help(h);
end
uiresume(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_eraser(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

n  = opt.radius;
if numel(n)==1, n = [n n n]; end

xi = opt.ijk(1)+(-n(1):n(1));
yi = opt.ijk(2)+(-n(2):n(2));
zi = opt.ijk(3)+(-n(3):n(3));

opt.mask(xi,yi,zi) = false;
opt.data = opt.data & opt.mask;

setappdata(h, 'opt', opt);
uiresume(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_erasercontour(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

ijk = opt.ijk;
n   = opt.radius*2+1;
if numel(n)==1, n = [n n n]; end

set(opt.handleseraser(1),'xdata',ijk(3)-n(3)/2+[0 0 n(3) n(3) 0]);
set(opt.handleseraser(1),'ydata',ijk(2)-n(2)/2+[0 n(2) n(2) 0 0]);
set(opt.handleseraser(2),'xdata',ijk(1)-n(1)/2+[0 0 n(1) n(1) 0]);
set(opt.handleseraser(2),'ydata',ijk(2)-n(2)/2+[0 n(2) n(2) 0 0]);
set(opt.handleseraser(3),'xdata',ijk(3)-n(3)/2+[0 0 n(3) n(3) 0]);
set(opt.handleseraser(3),'ydata',ijk(1)-n(1)/2+[0 n(1) n(1) 0 0]);

uiresume(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h, 'currentaxes');

xi  = opt.ijk(1);
yi  = opt.ijk(2);
zi  = opt.ijk(3);

dat1 = squeeze(opt.data(xi,:,:));
dat2 = opt.data(:,:,zi)';
dat3 = squeeze(opt.data(:,yi,:));

set(opt.handlesslice(1), 'CData', dat1);
set(opt.handlesslice(2), 'CData', dat2);
set(opt.handlesslice(3), 'CData', dat3);

crosshair([zi yi], 'handle', opt.handlescross(1,:));
crosshair([xi yi], 'handle', opt.handlescross(2,:));
crosshair([zi xi], 'handle', opt.handlescross(3,:));

cb_erasercontour(h);

set(h, 'currentaxes', curr_ax);

setappdata(h, 'opt', opt);
uiresume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonpress(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

cb_getposition(h);
opt = getappdata(h, 'opt');

seltype = get(h, 'selectiontype');
switch seltype
  case 'normal'
    % just update to new position, nothing else to be done here
    cb_redraw(h);
  case 'alt'
    cb_eraser(h);
    set(h, 'windowbuttonmotionfcn', @cb_tracemouse);
    opt = getappdata(h, 'opt');
    cb_redraw(h);
  otherwise
end

setappdata(h, 'opt', opt);

uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonrelease(h, eventdata)

seltype = get(h, 'selectiontype');
switch seltype
  case 'normal'
    % just update to new position, nothing else to be done here
  case 'alt'
    set(h, 'windowbuttonmotionfcn', '');  
  otherwise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_tracemouse(h, eventdata)

h   = getparent(h);
cb_getposition(h);
opt = getappdata(h, 'opt');

n  = opt.radius;
if numel(n)==1, n = [n n n]; end

xi = opt.ijk(1)+(-n(1):n(1)); xi(xi>opt.dim(1)) = []; xi(xi<1) = [];
yi = opt.ijk(2)+(-n(2):n(2)); yi(yi>opt.dim(2)) = []; yi(yi<1) = [];
zi = opt.ijk(3)+(-n(3):n(3)); zi(zi>opt.dim(3)) = []; zi(zi<1) = [];

opt.mask(xi,yi,zi) = 0;
opt.data = opt.data & opt.mask;

setappdata(h, 'opt', opt);
cb_redraw(h);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h,       'currentaxes');
pos     = get(curr_ax, 'currentpoint');

tag = get(curr_ax, 'tag');
switch tag
  case 'jk'
    opt.ijk([3,2]) = round(pos(1,1:2));
  case 'ji'
    opt.ijk([1,2]) = round(pos(1,1:2));
  case 'ik'
    opt.ijk([3,1]) = round(pos(1,1:2));
  otherwise
end
opt.ijk = min(opt.ijk, opt.dim);
opt.ijk = max(opt.ijk, [1 1 1]);

setappdata(h, 'opt', opt);
uiresume;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_cleanup(h, eventdata)

opt = getappdata(h, 'opt');
opt.quit = true;
setappdata(h, 'opt', opt);
uiresume
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = parseKeyboardEvent(eventdata)

key = eventdata.Key;

% handle possible numpad events (different for Windows and UNIX systems)
% NOTE: shift+numpad number does not work on UNIX, since the shift
% modifier is always sent for numpad events
if isunix()
  shiftInd = match_str(eventdata.Modifier, 'shift');
  if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
    % now we now it was a numpad keystroke (numeric character sent AND
    % shift modifier present)
    key = eventdata.Character;
    eventdata.Modifier(shiftInd) = []; % strip the shift modifier
  end
elseif ispc()
  if strfind(eventdata.Key, 'numpad')
    key = eventdata.Character;
  end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end
  
