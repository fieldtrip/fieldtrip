function data = ft_anonymizedata(cfg, data)

% FT_ANONYMIZEDATA clears the value of potentially identifying fields in
% the data and in the provenance information, i.e., it updates the data and
% the configuration structure and history that is maintained by FieldTrip
% in the cfg field.
%
% Use as
%   output = ft_anonymizedata(cfg, data)
% where data is any FieldTrip data structure and cfg is a configuration
% structure that should contain
%   cfg.keepnumeric = 'yes' or 'no', keep numeric fields (default = 'yes')
%   cfg.keepfield   = cell-array with strings, fields to keep (default = {})
%   cfg.removefield = cell-array with strings, fields to remove (default = {})
%   cfg.keepvalue   = cell-array with strings, values to keep (default = {})
%   cfg.removevalue = cell-array with strings, values to remove (default = {})
%
% The graphical user interface consists of a table that shows the name and
% value of each provenance element, and whether it should be kept or
% removed. Furthermore, it has a number of buttons:
%   - sort        specify which column is used for sorting
%   - apply       apply the current selection of 'keep' and 'remove' and hide the corresponding rows
%   - keep all    toggle all visibe rows to 'keep'
%   - remove all  toggle all visibe rows to 'keep'
%   - clear all   clear all visibe rows, i.e. neither 'keep' nor 'remove'
%   - quit        apply the current selection of 'keep' and 'remove' and exit
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile  = ...
%   cfg.outputfile  = ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_DEFACEVOLUME, FT_DEFACEMESH, FT_ANALYSISPIPELINE

% Copyright (C) 2014, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%  FieldTrip is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%  FieldTrip is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% get the options
cfg.keepfield   = ft_getopt(cfg, 'keepfield', {});
cfg.removefield = ft_getopt(cfg, 'removefield', {});
cfg.keepvalue   = ft_getopt(cfg, 'keepvalue', {});
cfg.removevalue = ft_getopt(cfg, 'removevalue', {});
cfg.keepnumeric = ft_getopt(cfg, 'keepnumeric', 'yes');

% determine the name and value of each element in the structure
[name, value] = splitstruct('data', data);

% we can rule out the numeric values as identifying
sel = cellfun(@ischar, value);
if istrue(cfg.keepnumeric)
  name  = name(sel);
  value = value(sel);
else
  % the numeric values are also to be judged by the end-user, but cannot be displayed easily
  % FIXME it would be possible to display single scalar values by converting them to a string
  value(~sel) = {'<numeric>'};
end

% do not bother with fields that are empty
sel = cellfun(@numel, value)>0;
name  = name(sel);
value = value(sel);

% all values are char, but some might be a char-array rather than a single string
sel = cellfun('size', value, 1)>1;
value(sel) = {'<multiline char array>'};

for i=1:length(value)
  % remove all non-printable characters
  sel = value{i}<32 | value{i}>126;
  value{i}(sel) = [];
end

keep = false(size(name));
for i=1:numel(cfg.keepfield)
  expression = sprintf('\\.%s$', cfg.keepfield{i});
  keep = keep | ~cellfun(@isempty, regexp(name, expression), 'uniformoutput', 1);
  expression = sprintf('\\.%s\\.', cfg.keepfield{i});
  keep = keep | ~cellfun(@isempty, regexp(name, expression), 'uniformoutput', 1);
  expression = sprintf('\\.%s\\(', cfg.keepfield{i});
  keep = keep | ~cellfun(@isempty, regexp(name, expression), 'uniformoutput', 1);
  expression = sprintf('\\.%s\\{', cfg.keepfield{i});
  keep = keep | ~cellfun(@isempty, regexp(name, expression), 'uniformoutput', 1);
end

keep = keep | ismember(value, cfg.keepvalue);

remove = false(size(name));
for i=1:numel(cfg.removefield)
  expression = sprintf('\\.%s$', cfg.removefield{i});
  remove = remove | ~cellfun(@isempty, regexp(name, expression), 'uniformoutput', 1);
  expression = sprintf('\\.%s\\.', cfg.removefield{i});
  remove = remove | ~cellfun(@isempty, regexp(name, expression), 'uniformoutput', 1);
  expression = sprintf('\\.%s\\(', cfg.removefield{i});
  remove = remove | ~cellfun(@isempty, regexp(name, expression), 'uniformoutput', 1);
  expression = sprintf('\\.%s\\{', cfg.removefield{i});
  remove = remove | ~cellfun(@isempty, regexp(name, expression), 'uniformoutput', 1);
end

remove = remove | ismember(value, cfg.removevalue);

% ensure there is no overlap
keep(remove) = false;

%% construct the graphical user interface

h = figure;
set(h, 'menuBar', 'none')
mp = get(0, 'MonitorPosition');
if size(mp,1)==1
  % there is only a single monitor, we can try to go fullscreen
  set(h, 'units', 'normalized', 'position', [0 0 1 1])
else
  set(h, 'units', 'normalized');
end

%% add the table to the GUI

t = uitable;
set(t, 'ColumnEditable', [true true false false]);
set(t, 'ColumnName', {'keep', 'remove', 'name', 'value'});
set(t, 'RowName', {});

%% add the buttons to the GUI

uicontrol('tag', 'button1', 'parent', h, 'units', 'pixels', 'style', 'popupmenu', 'string', {'remove', 'keep', 'name', 'value'}, 'userdata', 'sort', 'callback', @sort_cb);
uicontrol('tag', 'button2', 'parent', h, 'units', 'pixels', 'style', 'pushbutton', 'string', 'apply',       'userdata', 'a',  'callback', @keyboard_cb)
uicontrol('tag', 'button3', 'parent', h, 'units', 'pixels', 'style', 'pushbutton', 'string', 'keep all',    'userdata', 'ka', 'callback', @keyboard_cb)
uicontrol('tag', 'button4', 'parent', h, 'units', 'pixels', 'style', 'pushbutton', 'string', 'remove all',  'userdata', 'ra', 'callback', @keyboard_cb)
uicontrol('tag', 'button5', 'parent', h, 'units', 'pixels', 'style', 'pushbutton', 'string', 'clear all',   'userdata', 'ca', 'callback', @keyboard_cb)
uicontrol('tag', 'button6', 'parent', h, 'units', 'pixels', 'style', 'pushbutton', 'string', 'quit',        'userdata', 'q',  'callback', @keyboard_cb)

% use manual positioning of the buttons in pixel units
ft_uilayout(h, 'tag', 'button1', 'hpos', 20+(100+10)*0, 'vpos', 10, 'width', 100, 'height', 25);
ft_uilayout(h, 'tag', 'button2', 'hpos', 20+(100+10)*1, 'vpos', 10, 'width', 100, 'height', 25);
ft_uilayout(h, 'tag', 'button3', 'hpos', 20+(100+10)*2, 'vpos', 10, 'width', 100, 'height', 25);
ft_uilayout(h, 'tag', 'button4', 'hpos', 20+(100+10)*3, 'vpos', 10, 'width', 100, 'height', 25);
ft_uilayout(h, 'tag', 'button5', 'hpos', 20+(100+10)*4, 'vpos', 10, 'width', 100, 'height', 25);
ft_uilayout(h, 'tag', 'button6', 'hpos', 20+(100+10)*5, 'vpos', 10, 'width', 100, 'height', 25);

ft_uilayout(h, 'tag', 'button1', 'retag', 'buttongroup')
ft_uilayout(h, 'tag', 'button1', 'retag', 'buttongroup')
ft_uilayout(h, 'tag', 'button1', 'retag', 'buttongroup')
ft_uilayout(h, 'tag', 'button1', 'retag', 'buttongroup')
ft_uilayout(h, 'tag', 'button1', 'retag', 'buttongroup')
ft_uilayout(h, 'tag', 'button1', 'retag', 'buttongroup')

ft_uilayout(h, 'tag', 'buttongroup', 'BackgroundColor', [0.8 0.8 0.8]);

% this structure is passed around as appdata
info         = [];
info.table   = t;
info.name    = name;
info.value   = value;
info.keep    = keep;
info.remove  = remove;
info.hide    = false(size(name));
info.cleanup = false;
info.cfg     = cfg;

% this is consistent with the sort button
[dum, indx] = sort(~info.remove);
info.keep   = info.keep(indx);
info.remove = info.remove(indx);
info.name   = info.name(indx);
info.value  = info.value(indx);
info.hide   = info.hide(indx);

% these callbacks need the info appdata
setappdata(h, 'info', info);
set(h, 'CloseRequestFcn', 'delete(gcf)');
set(h, 'ResizeFcn', @resize_cb);

redraw_cb(h);
resize_cb(h);

while ~info.cleanup

  uiwait(h); % we only get part this point with abort or cleanup

  if ~ishandle(h)
    ft_error('aborted by user');
  end

  info = getappdata(h, 'info');

  if info.cleanup
    if ~all(xor(info.keep, info.remove))
      ft_warning('not all fields have been marked as "keep" or "remove"');
      info.cleanup = false;
    else
      delete(h);
    end
  end
end


fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
name = info.name(info.remove);
for i=1:length(name)
  str = sprintf('%s = ''removed by ft_anonymizedata'';', name{i});
  disp(str);
  eval(str);
end
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

% deal with the output
ft_postamble debug
ft_postamble previous data
ft_postamble provenance data
ft_postamble history data
ft_postamble savevar data

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end
end % function

function redraw_cb(h, eventdata)
h = getparent(h);
info = getappdata(h, 'info');
data = cat(2, num2cell(info.keep), num2cell(info.remove), info.name, info.value);
data = data(~info.hide,:);
set(info.table, 'data', data);
end % function

function keyboard_cb(h, eventdata)

if (isempty(eventdata) && ft_platform_supports('matlabversion',-Inf, '2014a')) || isa(eventdata, 'matlab.ui.eventdata.ActionData')
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

h = getparent(h);
info = getappdata(h, 'info');

data = get(info.table, 'data');

sel = info.keep & info.remove;
if any(sel)
  ft_warning('items that were marked both as "keep" and "remove" have been cleared');
  info.keep(sel) = false;
  info.remove(sel) = false;
end

info.keep  (~info.hide) = cell2mat(data(:,1));
info.remove(~info.hide) = cell2mat(data(:,2));

switch key
  case 'q'
    info.cleanup = true;
    setappdata(h, 'info', info); % store it immediately
    uiresume                     % resume from uiwait in the main function
  case 'a'
    info.hide(info.keep)   = true;
    info.hide(info.remove) = true;
  case 'ka'
    info.keep  (~info.hide) = true;
    info.remove(~info.hide) = false;
  case 'ra'
    info.keep  (~info.hide) = false;
    info.remove(~info.hide) = true;
  case 'ca'
    info.keep  (~info.hide) = false;
    info.remove(~info.hide) = false;
end
setappdata(h, 'info', info);
redraw_cb(h)
end % function

function sort_cb(h, eventdata)
h = getparent(h);
info = getappdata(h, 'info');
val = get(findobj(h, 'userdata', 'sort'), 'value');
str = get(findobj(h, 'userdata', 'sort'), 'string');
switch str{val}
  case 'remove'
    [dum, indx] = sort(~info.remove);
  case 'keep'
    [dum, indx] = sort(~info.keep);
  case 'name'
    [dum, indx] = sort(info.name);
  case 'value'
    [dum, indx] = sort(info.value);
end
info.keep   = info.keep(indx);
info.remove = info.remove(indx);
info.name   = info.name(indx);
info.value  = info.value(indx);
info.hide   = info.hide(indx);
setappdata(h, 'info', info);
redraw_cb(h);
end % function

function resize_cb(h, eventdata)
drawnow
h = getparent(h);
info = getappdata(h, 'info');
set(h, 'units', 'pixels');
siz = get(h, 'position');
% the 15 is for the vertical scrollbar on the right
set(info.table, 'units', 'normalized', 'position', [0.05 0.1 0.90 0.85]);
set(info.table, 'ColumnWidth', {50 50 300 600});
end
