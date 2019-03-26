function bnd = prepare_mesh_manual(cfg, mri)

% PREPARE_MESH_MANUAL is called by PREPARE_MESH and opens a GUI to manually
% select points/polygons in an mri dataset.
%
% It allows:
%   Visualization of 3d data in 3 different projections
%   Adjustment of brightness for every slice
%   Storage of the data points in an external .mat file
%   Retrieval of previously saved data points
%   Slice fast scrolling with keyboard arrows
%   Polygons or points selection/deselection
%
% See also PREPARE_MESH_SEGMENTATION, PREPARE_MESH_HEADSHAPE

% Copyrights (C) 2009, Cristiano Micheli & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
% $Id$

% FIXME: control slice's cmap referred to abs values
% FIXME: clean structure slicedata
% FIXME: check function assign3dpoints

global obj

mri = ft_checkdata(mri, 'datatype', {'volume', 'segmentation'}, 'hasunit', 'yes');

hasheadshape = isfield(cfg, 'headshape');
hasbnd       = isfield(cfg, 'bnd');  % FIXME why is this in cfg?
hasmri       = nargin>1;

% check the consistency of the input arguments
if hasheadshape && hasbnd
  ft_error('you should not specify cfg.headshape and cfg.bnd simultaneously');
end

% check the consistency of the input arguments
if ~hasmri
  % FIXME give a warning or so?
  mri.anatomy = [];
  mri.transform = eye(4);
else
  % ensure that it is double precision
  mri.anatomy = double(mri.anatomy);
end

% start with an empty boundary
bnd.pnt = [];
bnd.tri = [];

if hasheadshape
  if ~isempty(cfg.headshape)
    % start with the headshape
    [bnd.pos, bnd.tri] = headsurface([], [], 'headshape', cfg.headshape);
  end
elseif hasbnd
  if ~isempty(cfg.bnd)
    % start with the prespecified boundaries
    bnd = cfg.bnd;
  end
end

% creating the GUI
fig = figure;

% initialize some values
slicedata = []; points = bnd.pnt;
axs = floor(size(mri.anatomy,1)/2);

% initialize properties
set(fig, 'KeyPressFcn',@keypress);
set(fig, 'CloseRequestFcn', @cb_close);
setappdata(fig,'data',mri.anatomy);
setappdata(fig,'prop',{1,axs,0,[]});
setappdata(fig,'slicedata',slicedata);
setappdata(fig,'points',points);
setappdata(fig,'box',[]);

% add GUI elements
cb_creategui(gca);
cb_redraw(gca);
waitfor(fig);

% close sequence
global pnt
try
  tmp = pnt;
  clear global pnt
  pnt = tmp;
  clear tmp
  [tri] = projecttri(pnt);
  bnd.pnt = pnt;
  bnd.tri = tri;
catch
  bnd.pnt = [];
  bnd.tri = [];
end


function cb_redraw(hObject, eventdata, handles)
fig  = get(hObject, 'parent');
prop = getappdata(fig,'prop');
data = getappdata(fig,'data');
slicedata = getappdata(fig,'slicedata');
points = getappdata(fig,'points');

% draw image
try
  cla
  % i,j,k axes
  if (prop{1}==1) %jk
    imh=imagesc(squeeze(data(prop{2},:,:)));
    xlabel('k');
    ylabel('j');
    colormap bone
  elseif (prop{1}==2) %ik
    imh=imagesc(squeeze(data(:,prop{2},:)));
    xlabel('k');
    ylabel('i');
    colormap bone
  elseif (prop{1}==3) %ij
    imh=imagesc(squeeze(data(:,:,prop{2})));
    xlabel('j');
    ylabel('i');
    colormap bone
  end
  axis equal; axis tight
  brighten(prop{3});
  title([ 'slice ' num2str(prop{2}) ])
catch
end
if ~ishold,hold on,end
% draw points and lines
for zz = 1:length(slicedata)
  if(slicedata(zz).proj==prop{1})
    if(slicedata(zz).slice==prop{2})
      polygon = slicedata(zz).polygon;
      for kk=1:length(polygon)
        if ~isempty(polygon{kk})
          [xs,ys]=deal(polygon{kk}(:,1),polygon{kk}(:,2));
          plot(xs, ys, 'g.-');
        end
      end
    end
  end
end % end for

% draw other slices points
points = getappdata(fig,'points');
pnts   = round(points);
if ~isempty(pnts)
  % exclude polygons inside the same slice
  indx  = find(points(:,prop{1})==prop{2});
  tmp   = ones(1,size(points,1));
  tmp(indx)  = 0;
  pnts  = points(find(tmp),:);
  % calculate indexes of other proj falling in current slice
  indx2 = find(round(pnts(:,prop{1}))==prop{2});
  rest  = setdiff([1 2 3],prop{1});
  % plot the points
  for jj=1:length(indx2)
    [x,y] = deal(pnts(indx2(jj),rest(1)),pnts(indx2(jj),rest(2)));
    plot(y,x,'g*')
  end
end

% add blue cross markers in case of box selection
box        = getappdata(fig,'box');
point2mark = getappdata(fig,'point2mark');

if ~isempty(point2mark)
  plot(point2mark.x,point2mark.y,'marker','+')
end

if ~isempty(box)
  for kk=1:size(box,1)
    proj   = box(kk,1);
    slice  = box(kk,2);
    rowmin = box(kk,3);
    rowmax = box(kk,4);
    colmin = box(kk,5);
    colmax = box(kk,6);
    aaa = get(findobj(fig, 'color', 'g'));
    for ii=1:length(aaa)
      if ispolygon(aaa(ii))
        L = length(aaa(ii).YData);
        for jj=1:L
          cond = lt(aaa(ii).YData(jj),rowmax).*gt(aaa(ii).YData(jj),rowmin).* ...
            lt(aaa(ii).XData(jj),colmax).*gt(aaa(ii).XData(jj),colmin);
          if cond
            plot(aaa(ii).XData(jj),aaa(ii).YData(jj),'marker','+')
          end
        end
      elseif ispoint(aaa(ii))
        cond = lt(aaa(ii).YData,rowmax).*gt(aaa(ii).YData,rowmin).* ...
          lt(aaa(ii).XData,colmax).*gt(aaa(ii).XData,colmin);
        % the test after cond takes care the box doesnt propagate in 3D
        if (cond && (proj == prop{1}) && (slice == prop{2}))
          plot(aaa(ii).XData,aaa(ii).YData,'marker','+')
        end
      end
    end
  end
end

if get(findobj(fig, 'tag', 'toggle axes'), 'value')
  axis on
else
  axis off
end

function cb_creategui(hObject, eventdata, handles)

fig = get(hObject, 'parent');
% define the position of each GUI element
% constants
CONTROL_WIDTH  = 0.10;
CONTROL_HEIGHT = 0.075; 
CONTROL_HOFFSET = 0.8;
CONTROL_VOFFSET = 0.6;

% control buttons
uicontrol('tag', 'view', 'parent', fig, 'units', 'normalized', 'style', 'radiobutton', 'string', 'View', 'value', 1, 'callback', @which_task1);
ft_uilayout(fig, 'tag', 'view',  'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT, 'vpos', CONTROL_VOFFSET+4*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);

uicontrol('tag', 'ins', 'parent', fig, 'units', 'normalized', 'style', 'radiobutton', 'string', 'Ins', 'value', 0, 'callback', @which_task2);
ft_uilayout(fig, 'tag', 'ins',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT, 'vpos', CONTROL_VOFFSET+2.75*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'del', 'parent', fig, 'units', 'normalized', 'style', 'radiobutton', 'string', 'Del', 'value', 0, 'callback', @which_task3);
ft_uilayout(fig, 'tag', 'del',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT, 'vpos', CONTROL_VOFFSET+2.75*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+CONTROL_WIDTH);

uicontrol('tag', 'slice', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'slice', 'value', [], 'callback', []);
ft_uilayout(fig, 'tag', 'slice',  'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT, 'vpos', CONTROL_VOFFSET+1.25*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'min', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', '<', 'value', [], 'callback', @cb_btn);
ft_uilayout(fig, 'tag', 'min',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET+1*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'max', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', '>', 'value', [], 'callback', @cb_btn);
ft_uilayout(fig, 'tag', 'max',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET+1*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+CONTROL_WIDTH);

uicontrol('tag', 'brightness', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'brightness', 'value', [], 'callback', []);
ft_uilayout(fig, 'tag', 'brightness',  'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT, 'vpos', CONTROL_VOFFSET-0.75*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'plus', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', '+', 'value', [], 'callback', @cb_btn);
ft_uilayout(fig, 'tag', 'plus',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-1*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'minus', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', '-', 'value', [], 'callback', @cb_btn);
ft_uilayout(fig, 'tag', 'minus',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-1*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+CONTROL_WIDTH);

uicontrol('tag', 'toggle axes', 'parent', fig, 'units', 'normalized', 'style', 'checkbox', 'string', 'axes', 'value', 0, 'callback', @cb_redraw);
ft_uilayout(fig, 'tag', 'toggle axes',  'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT, 'vpos', CONTROL_VOFFSET-2.5*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);

uicontrol('tag', 'plane', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'plane', 'value', [], 'callback', []);
ft_uilayout(fig, 'tag', 'plane',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-3.5*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'one', 'parent', fig, 'units', 'normalized', 'style', 'radiobutton', 'string', 'jk', 'value', 0, 'callback', @cb_btnp1);
ft_uilayout(fig, 'tag', 'one',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-3.5*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+CONTROL_WIDTH);
uicontrol('tag', 'two', 'parent', fig, 'units', 'normalized', 'style', 'radiobutton', 'string', 'ik', 'value', 0, 'callback', @cb_btnp2);
ft_uilayout(fig, 'tag', 'two',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-4*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'three', 'parent', fig, 'units', 'normalized', 'style', 'radiobutton', 'string', 'ij', 'value', 0, 'callback', @cb_btnp3);
ft_uilayout(fig, 'tag', 'three',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-4*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+CONTROL_WIDTH);

uicontrol('tag', 'mesh', 'parent', fig, 'units', 'normalized', 'style', 'text', 'string', 'mesh', 'value', [], 'callback', []);
ft_uilayout(fig, 'tag', 'mesh',  'BackgroundColor', [0.8 0.8 0.8], 'width', 2*CONTROL_WIDTH, 'height', CONTROL_HEIGHT, 'vpos', CONTROL_VOFFSET-5.75*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);

uicontrol('tag', 'smoothm', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'smooth', 'value', [], 'callback', @smooth_mesh);
ft_uilayout(fig, 'tag', 'smoothm',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-6*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'viewm', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'view', 'value', [], 'callback', @view_mesh);
ft_uilayout(fig, 'tag', 'viewm',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-6*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+CONTROL_WIDTH);

uicontrol('tag', 'cancel', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'cancel', 'value', [], 'callback', @cancel_mesh);
ft_uilayout(fig, 'tag', 'cancel',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-6.5*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET);
uicontrol('tag', 'ok', 'parent', fig, 'units', 'normalized', 'style', 'pushbutton', 'string', 'ok', 'value', [], 'callback', @cb_btn);
ft_uilayout(fig, 'tag', 'ok',  'BackgroundColor', [0.8 0.8 0.8], 'width', CONTROL_WIDTH, 'height', CONTROL_HEIGHT/2, 'vpos', CONTROL_VOFFSET-6.5*CONTROL_HEIGHT, 'hpos', CONTROL_HOFFSET+CONTROL_WIDTH);


function cb_btn(hObject, eventdata, handles)
fig  = get(gca, 'parent');
prop = getappdata(fig,'prop');
slice = prop{2};
br    = prop{3};
datdim = size(getappdata(fig,'data'));
% p1 = [];p2 = [];p3 = [];
p1  = get(findobj('string','jk'),'value');
p2  = get(findobj('string','ik'),'value');
p3  = get(findobj('string','ij'),'value');
set(findobj('string','-'),'Value',prop{3});
set(findobj('string','+'),'Value',prop{3});
beta = get(findobj('string','-'),'Value');

if (p1)     %jk
  setappdata(fig,'prop',{1,round(datdim(1)/2),br,prop{4}});
  cb_redraw(gca);
elseif (p2) %ik
  setappdata(fig,'prop',{2,round(datdim(2)/2),br,prop{4}});
  cb_redraw(gca);
elseif (p3) %ij
  setappdata(fig,'prop',{3,round(datdim(3)/2),br,prop{4}});
  cb_redraw(gca);
end
if strcmp(get(hObject,'string'),'-')
  beta=beta-0.05;
  set(findobj('string','-'),'Value',beta);
  set(findobj('string','+'),'Value',beta);
  setappdata(fig,'prop',{prop{1},prop{2},beta,prop{4}});
  cb_redraw(gca);
elseif strcmp(get(hObject,'string'),'+')
  beta=beta+0.05;
  set(findobj('string','+'),'Value',beta);
  setappdata(fig,'prop',{prop{1},prop{2},beta,prop{4}});
  cb_redraw(gca);
end

if strcmp(get(hObject,'string'),'<')
  setappdata(fig,'prop',{prop{1},slice-1,0,prop{4}});
  cb_redraw(gca);
end
if strcmp(get(hObject,'string'),'>')
  setappdata(fig,'prop',{prop{1},slice+1,0,prop{4}});
  cb_redraw(gca);
end
if strcmp(get(hObject,'string'),'OK')
  close(fig)
end

function keypress(h, eventdata, handles, varargin)
fig   = get(gca, 'parent');
prop  = getappdata(fig,'prop');
dat   = guidata(gcbf);
key   = get(gcbf, 'CurrentCharacter');
slice = prop{2};
if  key
  switch key
    case 28
      setappdata(fig,'prop',{prop{1},slice-1,0,prop{4}});
      cb_redraw(gca);
    case 29
      setappdata(fig,'prop',{prop{1},slice+1,0,prop{4}});
      cb_redraw(gca);
    otherwise
  end
end
uiresume(h)

function build_polygon
fig   = get(gca, 'parent');
prop  = getappdata(fig,'prop');
data  = getappdata(fig,'data');
ss = size(data);
proj  = prop{1};
slice = prop{2};
thispolygon =1;
polygon{thispolygon} = zeros(0,2);
maskhelp = [ ...
  '------------------------------------------------------------------------\n' ...
  'specify polygons for masking the topographic interpolation\n' ...
  'press the right mouse button to add another point to the current polygon\n' ...
  'press backspace on the keyboard to remove the last point\n' ...
  'press "c" on the keyboard to close this polygon and start with another\n' ...
  'press "q" or ESC on the keyboard to continue\n' ...
  ];
again = 1;
if ~ishold,hold,end
fprintf(maskhelp);
fprintf('\n');

while again
  for i=1:length(polygon)
    fprintf('polygon %d has %d points\n', i, size(polygon{i},1));
  end
  
  [x, y, k] = ginput(1);
  okflag = 0;
  
  % check points do not fall out of image boundaries
  if get(findobj('string','Ins'),'value')
    if (proj == 1 && y<ss(2) && x<ss(3) && x>1 && y>1)
      okflag = 1;
    elseif (proj == 2 && y<ss(1) && x<ss(3) && x>1 && y>1)
      okflag = 1;
    elseif (proj == 3 && y<ss(1) && x<ss(2) && x>1 && y>1)
      okflag = 1;
    else
      okflag = 0;
    end
  end
  
  if okflag
    switch lower(k)
      case 1
        polygon{thispolygon} = cat(1, polygon{thispolygon}, [x y]);
        % add the last line segment to the figure
        if size(polygon{thispolygon},1)>1
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
        end
        plot(x, y, 'g.-');
        
      case 8 % backspace
        if size(polygon{thispolygon},1)>0
          % redraw existing polygons
          cb_redraw(gca);
          % remove the last point of current polygon
          polygon{thispolygon} = polygon{thispolygon}(1:end-1,:);
          for i=1:length(polygon)
            x = polygon{i}(:,1);
            y = polygon{i}(:,2);
            if i~=thispolygon
              % close the polygon in the figure
              x(end) = x(1);
              y(end) = y(1);
            end
            set(gca,'nextplot','new')
            plot(x, y, 'g.-');
          end
        end
        
      case 'c'
        if size(polygon{thispolygon},1)>0
          % close the polygon
          polygon{thispolygon}(end+1,:) = polygon{thispolygon}(1,:);
          % close the polygon in the figure
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
          plot(x, y, 'g.-');
          % switch to the next polygon
          thispolygon = thispolygon + 1;
          polygon{thispolygon} = zeros(0,2);
          slicedata = getappdata(fig,'slicedata');
          m = length(slicedata);
          slicedata(m+1).proj    = prop{1};
          slicedata(m+1).slice   = prop{2};
          slicedata(m+1).polygon = polygon;
          setappdata(fig,'slicedata',slicedata);
        end
        
      case {'q', 27}
        if size(polygon{thispolygon},1)>0
          % close the polygon
          polygon{thispolygon}(end+1,:) = polygon{thispolygon}(1,:);
          % close the polygon in the figure
          x = polygon{i}([end-1 end],1);
          y = polygon{i}([end-1 end],2);
          plot(x, y, 'g.-');
        end
        again = 0;
        %set(gcf,'WindowButtonDownFcn','');
        slicedata = getappdata(fig,'slicedata');
        m = length(slicedata);
        if ~isempty(polygon{thispolygon})
          slicedata(m+1).proj    = prop{1};
          slicedata(m+1).slice   = prop{2};
          slicedata(m+1).polygon = polygon;
          setappdata(fig,'slicedata',slicedata);
        end
      otherwise
        ft_warning('invalid button (%d)', k);
    end
    
  end
  assign_points3d
end

function cb_btn_dwn(hObject, eventdata, handles)
build_polygon
uiresume

function cb_btn_dwn2(hObject, eventdata, handles)
global obj
erase_points2d(obj)
uiresume


function erase_points2d (obj)
fig       = get(gca, 'parent');
prop      = getappdata(fig,'prop');
slicedata = getappdata(fig,'slicedata');
pntsslice = [];

% 2d slice selected box boundaries
if strcmp(get(obj,'string'),'box')
  [x, y] = ft_select_box;
  rowmin = min(y);    rowmax = max(y);
  colmin = min(x);    colmax = max(x);
  box = getappdata(fig,'box');
  box = [box; prop{1} prop{2} rowmin rowmax colmin colmax];
  setappdata(fig,'box',box);
else
  % converts lines to points
  pos = slicedata2pnts(prop{1},prop{2});
  tmp = ft_select_point(pos);
  point2mark.x = tmp(1); point2mark.y = tmp(2);
  setappdata(fig,'point2mark',point2mark);
end

maskhelp = [ ...
  '------------------------------------------------------------------------\n' ...
  'Delete points:\n' ...
  'press "c" on the keyboard to undo all selections\n' ...
  'press "d" or DEL on the keyboard to delete all selections\n' ...
  ];
fprintf(maskhelp);
fprintf('\n');
again = 1;
if ~ishold,hold,end
cb_redraw(gca);

% waits for a key press
while again
  [junk1, junk2, k] = ginput(1);
  
  switch lower(k)
    case 'c'
      if strcmp(get(obj,'string'),'box')
        % undoes all the selections, but only in the current slice
        ind = [];
        for ii=1:size(box,1)
          if box(ii,1)==prop{1} && box(ii,2)==prop{2}
            ind = [ind,ii];
          end
        end
        box(ind,:) = [];
        setappdata(fig,'box',box);
        cb_redraw(gca);
        again = 0;
        set(gcf,'WindowButtonDownFcn','');
      else
        setappdata(fig,'point2mark',[]);
        cb_redraw(gca);
        again = 0;
        set(gcf,'WindowButtonDownFcn','');
      end
    case {'d', 127}
      if strcmp(get(obj,'string'),'box')
        update_slicedata
        % deletes only the points selected in the current slice
        ind = [];
        for ii=1:size(box,1)
          if box(ii,1)==prop{1} && box(ii,2)==prop{2}
            ind = [ind,ii];
          end
        end
        box(ind,:) = [];
        setappdata(fig,'box',box);
        cb_redraw(gca);
        again = 0;
        set(gcf,'WindowButtonDownFcn','');
      else
        update_slicedata
        setappdata(fig,'point2mark',[]);
        cb_redraw(gca);
        again = 0;
        set(gcf,'WindowButtonDownFcn','');
      end
    otherwise
  end
end
assign_points3d

function assign_points3d
fig       = get(gca, 'parent');
slicedata = getappdata(fig,'slicedata');

points = [];
for zz = 1:length(slicedata)
  slice   = slicedata(zz).slice;
  proj    = slicedata(zz).proj;
  polygon = slicedata(zz).polygon;
  for kk=1:length(polygon)
    if ~isempty(polygon{kk})
      [xs,ys]=deal(polygon{kk}(:,1),polygon{kk}(:,2));
      if (proj==1)     %jk
        tmp    = [ slice*ones(length(xs),1),ys(:),xs(:) ];
        points = [points; unique(tmp,'rows')];
      elseif (proj==2) %ik
        tmp    = [ ys(:),slice*ones(length(xs),1),xs(:) ];
        points = [points; unique(tmp,'rows')];
      elseif (proj==3) %ij
        tmp    = [ ys(:),xs(:),slice*ones(length(xs),1) ];
        points = [points; unique(tmp,'rows')];
      end
    end
  end
end
setappdata(fig,'points',points);

function update_slicedata
fig       = get(gca, 'parent');
prop      = getappdata(fig,'prop');
slicedata = getappdata(fig,'slicedata');
box       = getappdata(fig,'box');
point2mark = getappdata(fig,'point2mark');
% case of box selection
if ~isempty(box) && ~isempty(slicedata)
  for zz = 1:length(slicedata)
    slice   = slicedata(zz).slice;
    proj    = slicedata(zz).proj;
    polygon = slicedata(zz).polygon;
    for uu=1:size(box,1)
      if (box(uu,1)==proj) && (box(uu,2)==slice)
        for kk=1:length(polygon)
          if ~isempty(polygon{kk})
            [xs,ys] = deal(polygon{kk}(:,1),polygon{kk}(:,2));
            tmpind = lt(ys,ones(length(xs),1)*box(uu,4)).*gt(ys,ones(length(xs),1)*box(uu,3)).* ...
              lt(xs,ones(length(xs),1)*box(uu,6)).*gt(xs,ones(length(xs),1)*box(uu,5));
            slicedata(zz).polygon{kk} = [ xs(find(~tmpind)),ys(find(~tmpind)) ];
          end
        end
      end
    end
  end
  % case of point selection
else
  if ~isempty(slicedata)
    for zz = 1:length(slicedata)
      slice   = slicedata(zz).slice;
      proj    = slicedata(zz).proj;
      polygon = slicedata(zz).polygon;
      for kk=1:length(polygon)
        if ~isempty(polygon{kk})
          [xs,ys] = deal(polygon{kk}(:,1),polygon{kk}(:,2));
          tmpind = eq(xs,point2mark.x).*eq(ys,point2mark.y);
          slicedata(zz).polygon{kk} = [ xs(find(~tmpind)),ys(find(~tmpind)) ];
        end
      end
    end
  end
end
setappdata(fig,'slicedata',slicedata)

% FIXME: obsolete: replaced by headshape at the beginning
function smooth_mesh(hObject, eventdata, handles)
fig = get(gca,'parent');
disp('not yet implemented')

% FIXME: make it work for pre-loaded meshes
function view_mesh(hObject, eventdata, handles)
fig  = get(gca, 'parent');
assign_points3d
pnt  = getappdata(fig,'points');
pnt_ = pnt;
if ~isempty(pnt)
  tri = projecttri(pnt);
  bnd.pnt = pnt_;
  bnd.tri = tri;
  slicedata = getappdata(fig,'slicedata');
  figure
  ft_plot_mesh(bnd,'vertexcolor','k');
end

function cancel_mesh(hObject, eventdata, handles)
fig  = get(gca, 'parent');
close(fig)

function cb_close(hObject, eventdata, handles)
% get the points from the figure
fig    = hObject;
global pnt
pnt   = getappdata(fig, 'points');
set(fig, 'CloseRequestFcn', @delete);
delete(fig);

function cb_btnp1(hObject, eventdata, handles)
set(findobj('string','jk'),'value',1);
set(findobj('string','ik'),'value',0);
set(findobj('string','ij'),'value',0);
cb_btn(hObject);

function cb_btnp2(hObject, eventdata, handles)
set(findobj('string','jk'),'value',0);
set(findobj('string','ik'),'value',1);
set(findobj('string','ij'),'value',0);
cb_btn(hObject);

function cb_btnp3(hObject, eventdata, handles)
set(findobj('string','jk'),'value',0);
set(findobj('string','ik'),'value',0);
set(findobj('string','ij'),'value',1);
cb_btn(hObject);

function which_task1(hObject, eventdata, handles)
set(gcf,'WindowButtonDownFcn','');
set(findobj('string','View'),'value',1);
set(findobj('string','Ins'), 'value',0);
set(findobj('string','Del'), 'value',0);
set(gcf,'WindowButtonDownFcn','');
change_panel

function which_task2(hObject, eventdata, handles)
set(findobj('string','View'),'value',0);
set(findobj('string','Ins'), 'value',1);
set(findobj('string','Del'), 'value',0);
set(gcf,'WindowButtonDownFcn',@cb_btn_dwn);
change_panel

function which_task3(hObject, eventdata, handles)
set(gcf,'WindowButtonDownFcn','');
set(findobj('string','View'),'value',0);
set(findobj('string','Ins'), 'value',0);
set(findobj('string','Del'), 'value',1);
change_panel

function change_panel
% affect panel visualization
w1 = get(findobj('string','View'),'value');
w2 = get(findobj('string','Ins'), 'value');
w3 = get(findobj('string','Del'), 'value');

if w1
  set(findobj('tag','slice'),'string','slice');
  set(findobj('tag','min'),'string','<');
  set(findobj('tag','max'),'string','>');
  set(findobj('tag','min'),'style','pushbutton');
  set(findobj('tag','max'),'style','pushbutton');
  set(findobj('tag','min'),'callback',@cb_btn);
  set(findobj('tag','max'),'callback',@cb_btn);
  set(findobj('tag','mesh'),'visible','on');
  set(findobj('tag','openm'),'visible','on');
  set(findobj('tag','viewm'),'visible','on');
  set(findobj('tag','brightness'),'visible','on');
  set(findobj('tag','plus'),'visible','on');
  set(findobj('tag','minus'),'visible','on');
  set(findobj('tag','toggle axes'),'visible','on');
  set(findobj('tag','plane'),'visible','on');
  set(findobj('tag','one'),'visible','on');
  set(findobj('tag','two'),'visible','on');
  set(findobj('tag','three'),'visible','on');
elseif w2
  set(findobj('tag','slice'),'string','select points');
  set(findobj('tag','min'),'string','close');
  set(findobj('tag','max'),'string','next');
  set(findobj('tag','min'),'style','pushbutton');
  set(findobj('tag','max'),'style','pushbutton');
  set(findobj('tag','min'),'callback',@tins_close);
  set(findobj('tag','max'),'callback',@tins_next);
  set(findobj('tag','mesh'),'visible','off');
  set(findobj('tag','openm'),'visible','off');
  set(findobj('tag','viewm'),'visible','off');
  set(findobj('tag','brightness'),'visible','off');
  set(findobj('tag','plus'),'visible','off');
  set(findobj('tag','minus'),'visible','off');
  set(findobj('tag','toggle axes'),'visible','off');
  set(findobj('tag','plane'),'visible','off');
  set(findobj('tag','one'),'visible','off');
  set(findobj('tag','two'),'visible','off');
  set(findobj('tag','three'),'visible','off');
elseif w3
  set(findobj('tag','slice'),'string','select points');
  set(findobj('tag','min'),'string','box');
  set(findobj('tag','max'),'string','point');
  set(findobj('tag','min'),'style','pushbutton');
  set(findobj('tag','max'),'style','pushbutton');
  set(findobj('tag','min'),'callback',@tdel_box);
  set(findobj('tag','max'),'callback',@tdel_single);
  set(findobj('tag','mesh'),'visible','off');
  set(findobj('tag','openm'),'visible','off');
  set(findobj('tag','viewm'),'visible','off');
  set(findobj('tag','brightness'),'visible','off');
  set(findobj('tag','plus'),'visible','off');
  set(findobj('tag','minus'),'visible','off');
  set(findobj('tag','toggle axes'),'visible','off');
  set(findobj('tag','plane'),'visible','off');
  set(findobj('tag','one'),'visible','off');
  set(findobj('tag','two'),'visible','off');
  set(findobj('tag','three'),'visible','off');
end

function tins_close(hObject, eventdata, handles)
set(gcf,'WindowButtonDownFcn','');

function tins_next(hObject, eventdata, handles)
set(gcf,'WindowButtonDownFcn','');


function tdel_box(hObject, eventdata, handles)
global obj
obj = gco;
set(gcf,'WindowButtonDownFcn',@cb_btn_dwn2);

function tdel_single(hObject, eventdata, handles)
global obj
obj = gco;
set(gcf,'WindowButtonDownFcn',@cb_btn_dwn2);

% accessory functions:
function [pnts] = slicedata2pnts(proj,slice)
fig       = get(gca, 'parent');
slicedata = getappdata(fig,'slicedata');
pnts = [];
for i=1:length(slicedata)
  if slicedata(i).slice==slice && slicedata(i).proj == proj
    for j=1:length(slicedata(i).polygon)
      if ~isempty(slicedata(i).polygon{j})
        tmp = slicedata(i).polygon{j};
        xs = tmp(:,1); ys = tmp(:,2);
        pnts = [pnts;[unique([xs,ys],'rows')]];
      end
    end
  end
end

function h = ispolygon(x)
h = length(x.XData)>1;

function h = ispoint(x)
h = length(x.XData)==1;

function [T,P] = points2param(pnt)
x = pnt(:,1); y = pnt(:,2); z = pnt(:,3);
[phi,theta,R] = cart2sph(x,y,z);
theta = theta + pi/2;
phi   = phi + pi;
T     = theta(:);
P     = phi(:);

function [bnd2,a,Or] = spherical_harmonic_mesh(bnd, nl)

% SPHERICAL_HARMONIC_MESH realizes the smoothed version of a mesh contained in the first
% argument bnd. The boundary argument (bnd) contains typically 2 fields
% called .pnt and .tri referring to vertices and triangulation of a mesh.
% The degree of smoothing is given by the order in the second input: nl
%
% Use as
%   bnd2 = spherical_harmonic_mesh(bnd, nl);

% nl    = 12; % from van 't Ent paper

bnd2 = [];

% calculate midpoint
Or = mean(bnd.pnt);
% rescale all points
x = bnd.pnt(:,1) - Or(1);
y = bnd.pnt(:,2) - Or(2);
z = bnd.pnt(:,3) - Or(3);
X = [x(:), y(:), z(:)];

% convert points to parameters
[T,P] = points2param(X);

% basis function
B     = shlib_B(nl, T, P);
Y     = B'*X;

% solve the linear system
a     = pinv(B'*B)*Y;

% build the surface
Xs = zeros(size(X));
for l = 0:nl-1
  for m = -l:l
    if m<0
      Yml = shlib_Yml(l, abs(m), T, P);
      Yml = (-1)^m*conj(Yml);
    else
      Yml = shlib_Yml(l, m, T, P);
    end
    indx = l^2 + l + m+1;
    Xs = Xs + Yml*a(indx,:);
  end
  fprintf('%d of %d\n', l, nl);
end
% Just take the real part
Xs = real(Xs);
% reconstruct the matrices for plotting.
xs = reshape(Xs(:,1), size(x,1), size(x,2));
ys = reshape(Xs(:,2), size(x,1), size(x,2));
zs = reshape(Xs(:,3), size(x,1), size(x,2));

bnd2.pnt = [xs(:)+Or(1) ys(:)+Or(2) zs(:)+Or(3)];
[bnd2.tri] = projecttri(bnd.pnt);

function B = shlib_B(nl, theta, phi)
% function B = shlib_B(nl, theta, phi)
%
% Constructs the matrix of basis functions
% where b(i,l^2 + l + m) = Yml(theta, phi)
%
% See also: shlib_Yml.m, shlib_decomp.m, shlib_gen_shfnc.m
%
% Dr. A. I. Hanna (2006).
B = zeros(length(theta), nl^2+2*nl+1);
for l = 0:nl-1
  for m = -l:l
    if m<0
      Yml = shlib_Yml(l, abs(m), theta, phi);
      Yml = (-1)^m*conj(Yml);
    else
      Yml = shlib_Yml(l, m, theta, phi);
    end
    indx = l^2 + l + m+1;
    B(:, indx) = Yml;
  end
end
return;

function Yml = shlib_Yml(l, m, theta, phi)
% function Yml = shlib_Yml(l, m, theta, phi)
%
% MATLAB function that takes a given order and degree, and the matrix of
% theta and phi and constructs a spherical harmonic from these. The
% analogue in the 1D case would be to give a particular frequency.
%
% Inputs:
%  m - order of the spherical harmonic
%  l - degree of the spherical harmonic
%  theta - matrix of polar coordinates \theta \in [0, \pi]
%  phi - matrix of azimuthal cooridinates \phi \in [0, 2\pi)
%
% Example:
%
% [x, y] = meshgrid(-1:.1:1);
% z = x.^2 + y.^2;
% [phi,theta,R] = cart2sph(x,y,z);
% Yml = spharm_aih(2,2, theta(:), phi(:));
%
% See also: shlib_B.m, shlib_decomp.m, shlib_gen_shfnc.m
%
% Dr. A. I. Hanna (2006)
Pml=legendre(l,cos(theta));
if l~=0
  Pml=squeeze(Pml(m+1,:,:));
end
Pml = Pml(:);
% Yml = sqrt(((2*l+1)/4).*(factorial(l-m)/factorial(l+m))).*Pml.*exp(sqrt(-1).*m.*phi);
% new:
Yml = sqrt(((2*l+1)/(4*pi)).*(factorial(l-m)/factorial(l+m))).*Pml.*exp(sqrt(-1).*m.*phi);
return;
