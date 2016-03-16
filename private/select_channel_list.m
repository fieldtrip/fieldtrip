function [select] = select_channel_list(label, select, titlestr)

% SELECT_CHANNEL_LIST presents a dialog for selecting multiple elements from a cell
% array with strings, such as the labels of EEG channels. The dialog presents two
% columns with the mechanism to move elements between the two.
%
% Use as
%   final = select_channel_list(label, initial, titlestr)
% or 
%   final = select_channel_list(layout, initial, titlestr)
%
% with
%   initial   = indices of channels that are initially selected
%   label     = cell array with channel labels
%   layout    = structure with 2D channel layout, see FT_PREPARE_LAYOUT
%   titlestr  = title for dialog (optional, default is 'Select')
% and
%   final     = indices of the selected channels
%
% If the user presses cancel or closes the figure, the initial selection will be
% returned.

% Copyright (C) 2003-2016, Robert Oostenveld
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

if nargin<3
  % it can be channels, but also other items
  titlestr = 'Select';
end

% ensure that it is a row array
select = select(:)';

if isstruct(label) && all(isfield(label, {'label', 'pos', 'width', 'height'}))
  % implement the selection GUI using the topography of the channels
  layout = label; clear label;
  
  % the callback sequence is SELECT_CHANNEL_LIST(this) -> FT_SELECT_CHANNEL-> FT_SELECT_RANGE -> LAYOUT_SELECT(subfunction)
  
  pos      = get(0, 'DefaultFigurePosition');
  pos(3:4) = [600 400];
  
  dlg = figure('WindowStyle', 'modal', 'NumberTitle', 'off', 'Name', titlestr, 'Position', pos);
  set(dlg, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears
  
  userdata.layout   = layout;
  userdata.h1       = [];
  userdata.h2       = [];
  userdata.select   = select;
  userdata.unselect = setdiff(1:length(layout.label), select);
  set(dlg, 'userdata', userdata);
  
  uicontrol(dlg, 'style', 'pushbutton', 'position', [020  10    80  20], 'string', 'select all',   'callback', @layout_select_all);
  uicontrol(dlg, 'style', 'pushbutton', 'position', [120  10    80  20], 'string', 'select none',  'callback', @layout_select_none);
  uicontrol(dlg, 'style', 'pushbutton', 'position', [400  10    80  20], 'string', 'Cancel',       'callback', 'close');
  uicontrol(dlg, 'style', 'pushbutton', 'position', [500  10    80  20], 'string', 'OK',           'callback', 'uiresume');
  layout_redraw(dlg);
  
  set(dlg, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', {@layout_select, dlg}, 'event', 'WindowButtonUpFcn'});
  set(dlg, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', {@layout_select, dlg}, 'event', 'WindowButtonDownFcn'});
  set(dlg, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', {@layout_select, dlg}, 'event', 'WindowButtonMotionFcn'});
  
else
  % implement the selection GUI as a dialog with two columns
  
  pos      = get(0, 'DefaultFigurePosition');
  pos(3:4) = [290 300];
  dlg      = dialog('Name', titlestr, 'Position', pos);
  set(gca, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears
  
  userdata.label    = label;
  userdata.select   = select;
  userdata.unselect = setdiff(1:length(label), select);
  set(dlg, 'userdata', userdata);
  
  uicontrol(dlg, 'style', 'text',       'position', [ 10 240+20 80  20], 'string', 'unselected');
  uicontrol(dlg, 'style', 'text',       'position', [200 240+20 80  20], 'string', 'selected  ');
  uicontrol(dlg, 'style', 'listbox',    'position', [ 10  40+20 80 200], 'min', 0, 'max', 2, 'tag', 'lbunsel')
  uicontrol(dlg, 'style', 'listbox',    'position', [200  40+20 80 200], 'min', 0, 'max', 2, 'tag', 'lbsel')
  uicontrol(dlg, 'style', 'pushbutton', 'position', [105 175+20 80  20], 'string', 'add all >'   , 'callback', @label_addall);
  uicontrol(dlg, 'style', 'pushbutton', 'position', [105 145+20 80  20], 'string', 'add >'       , 'callback', @label_add);
  uicontrol(dlg, 'style', 'pushbutton', 'position', [105 115+20 80  20], 'string', '< remove'    , 'callback', @label_remove);
  uicontrol(dlg, 'style', 'pushbutton', 'position', [105  85+20 80  20], 'string', '< remove all', 'callback', @label_removeall);
  uicontrol(dlg, 'style', 'pushbutton', 'position', [ 55  10    80  20], 'string', 'Cancel',       'callback', 'close');
  uicontrol(dlg, 'style', 'pushbutton', 'position', [155  10    80  20], 'string', 'OK',           'callback', 'uiresume');
  label_redraw(dlg);
end % if layout

% wait untill the dialog is closed or the user presses OK/Cancel
uiwait(dlg);
if ishandle(dlg)
  % the user pressed OK, return the selection from the dialog
  userdata = get(dlg, 'userdata');
  select = userdata.select;
  close(dlg);
  return
else
  % the user pressed Cancel or closed the dialog, return the initial selection
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_redraw(h)
userdata = get(h, 'userdata');
set(findobj(h, 'tag', 'lbsel'  ), 'string', userdata.label(userdata.select));
set(findobj(h, 'tag', 'lbunsel'), 'string', userdata.label(userdata.unselect));
% set the active element in the select listbox, based on the previous active element
tmp = min(get(findobj(h, 'tag', 'lbsel'), 'value'));
tmp = min(tmp, length(get(findobj(h, 'tag', 'lbsel'), 'string')));
if isempty(tmp) | tmp==0
  tmp = 1;
end
set(findobj(h, 'tag', 'lbsel'  ), 'value', tmp);
% set the active element in the unselect listbox, based on the previous active element
tmp = min(get(findobj(h, 'tag', 'lbunsel'), 'value'));
tmp = min(tmp, length(get(findobj(h, 'tag', 'lbunsel'), 'string')));
if isempty(tmp) || tmp==0
  tmp = 1;
end
set(findobj(h, 'tag', 'lbunsel'  ), 'value', tmp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_addall(h, eventdata, handles, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
userdata.select   = 1:length(userdata.label);
userdata.unselect = [];
set(findobj(h, 'tag', 'lbunsel'  ), 'value', 1);
set(h, 'userdata', userdata);
label_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_removeall(h, eventdata, handles, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
userdata.unselect = 1:length(userdata.label);
userdata.select   = [];
set(findobj(h, 'tag', 'lbsel'  ), 'value', 1);
set(h, 'userdata', userdata);
label_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_add(h, eventdata, handles, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
if ~isempty(userdata.unselect)
  add = userdata.unselect(get(findobj(h, 'tag', 'lbunsel'  ), 'value'));
  userdata.select   = sort([userdata.select add]);
  userdata.unselect = sort(setdiff(userdata.unselect, add));
  set(h, 'userdata', userdata);
  label_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label_remove(h, eventdata, handles, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
if ~isempty(userdata.select)
  remove = userdata.select(get(findobj(h, 'tag', 'lbsel'  ), 'value'));
  userdata.select   = sort(setdiff(userdata.select, remove));
  userdata.unselect = sort([userdata.unselect remove]);
  set(h, 'userdata', userdata);
  label_redraw(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout_select_all(h, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
userdata.select   = 1:numel(userdata.layout.label);
userdata.unselect = [];
set(h, 'userdata', userdata);
layout_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout_select_none(h, varargin)
h = get(h, 'parent');
userdata = get(h, 'userdata');
userdata.select   = [];
userdata.unselect = 1:numel(userdata.layout.label);
set(h, 'userdata', userdata);
layout_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout_select(list, h)
% this is a callback executed by FT_SELECT_CHANNEL and FT_SELECT_RANGE
userdata = get(h, 'userdata');
userdata.select   = find( ismember(userdata.layout.label, list));
userdata.unselect = find(~ismember(userdata.layout.label, list));
set(h, 'userdata', userdata);
layout_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function layout_redraw(h)
userdata = get(h, 'userdata');
layout = userdata.layout;
% don't clear the whole figure, since it contains some buttons
% don't clear the whole axes, since it may contain multiple selected ranges
% only clear the plotted layout
delete(userdata.h1);
delete(userdata.h2);
% plot the selected and unselected channels differently
layout1 = [];
layout1.pos    = layout.pos(userdata.select,:);
layout1.width  = layout.width(userdata.select);
layout1.height = layout.height(userdata.select);
layout1.label  = layout.label(userdata.select);
userdata.h1 = ft_plot_lay(layout1, 'point', true, 'label', true, 'box', false, 'figure', h, 'labelcolor', [0 0 0]);
layout2 = [];
layout2.pos    = layout.pos(userdata.unselect,:);
layout2.width  = layout.width(userdata.unselect);
layout2.height = layout.height(userdata.unselect);
layout2.label  = layout.label(userdata.unselect);
userdata.h2 = ft_plot_lay(layout2, 'point', true, 'label', true, 'box', false, 'figure', h, 'labelcolor', [1 1 1]/2);
axis([min(layout.pos(:,1)) max(layout.pos(:,1)) min(layout.pos(:,2)) max(layout.pos(:,2))])
set(h, 'userdata', userdata);
% add the channel information to the figure, this is used in the callbacks
info          = guidata(h);
info.x        = layout.pos(:, 1);
info.y        = layout.pos(:, 2);
info.label    = layout.label;
guidata(h, info);

