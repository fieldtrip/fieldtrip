function [select] = select_channel_list(label, select, titlestr)

% SELECT_CHANNEL_LIST presents a dialog for selecting multiple elements
% from a cell-array with strings, such as the labels of EEG channels.
% The dialog presents two columns with an add and remove mechanism.
% 
% select = select_channel_list(label, initial, titlestr)
% 
% with 
%   initial indices of channels that are initially selected 
%   label   cell-array with channel labels (strings)
%   titlestr    title for dialog (optional)
% and
%   select  indices of selected channels
%
% If the user presses cancel, the initial selection will be returned.

% Copyright (C) 2003, Robert Oostenveld
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
  titlestr = 'Select';
end

pos      = get(0,'DefaultFigurePosition');
pos(3:4) = [290 300];
dlg      = dialog('Name', titlestr, 'Position', pos);
drawnow
set(dlg, 'Visible', 'off'); % explicitly turn the axis off, as it sometimes appears

select            = select(:)';     % ensure that it is a row array
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
if isempty(tmp) || tmp==0
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
function label_remove(h, eventdata, handles, varargin); 
h = get(h, 'parent');
userdata = get(h, 'userdata');
if ~isempty(userdata.select)
  remove = userdata.select(get(findobj(h, 'tag', 'lbsel'  ), 'value'));
  userdata.select   = sort(setdiff(userdata.select, remove));
  userdata.unselect = sort([userdata.unselect remove]);
  set(h, 'userdata', userdata);
  label_redraw(h);
end

