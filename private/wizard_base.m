function h = wizard_gui(filename)

% This is the low level wizard function. It evaluates the MATLAB content
% in the workspace of the calling function. To prevent overwriting
% variables in the BASE workspace, this function should be called from a
% wrapper function. The wrapper function whoudl pause execution untill the
% wizard figure is deleted.

% Copyright (C) 2007, Robert Oostenveld
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

% create a new figure
h = figure('Name','Wizard',...
  'NumberTitle','off',...
  'MenuBar','none',...
  'KeyPressFcn', @cb_keyboard,...
  'resizeFcn', @cb_resize,...
  'closeRequestFcn', @cb_cancel);
% movegui(h,'center');
% set(h, 'ToolBar', 'figure');

if nargin>0
  filename = which(filename);
  [p, f, x] = fileparts(filename);
  data          = [];
  data.path     = p;
  data.file     = f;
  data.ext      = x;
  data.current  = 0;
  data.script   = script_parse(filename);
  guidata(h, data);   % attach the data to the GUI figure
else
  data          = [];
  data.path     = [];
  data.file     = [];
  data.ext      = [];
  data.current  = 0;
  data.script   = script_parse([]);
  guidata(h, data);   % attach the data to the GUI figure
  cb_load(h);
end
cb_show(h);           % show the GUI figure
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)
h = parentfig(h);
dbstack
if     isequal(eventdata.Key, 'o') && isequal(eventdata.Modifier, {'control'})
  cb_load(h);
elseif isequal(eventdata.Key, 's') && isequal(eventdata.Modifier, {'control'})
  cb_save(h);
elseif isequal(eventdata.Key, 'e') && isequal(eventdata.Modifier, {'control'})
  cb_edit(h);
elseif isequal(eventdata.Key, 'p') && isequal(eventdata.Modifier, {'control'})
  cb_prev(h);
elseif isequal(eventdata.Key, 'n') && isequal(eventdata.Modifier, {'control'})
  cb_next(h);
elseif isequal(eventdata.Key, 'q') && isequal(eventdata.Modifier, {'control'})
  cb_cancel(h);
elseif isequal(eventdata.Key, 'x') && isequal(eventdata.Modifier, {'control'})
  % FIXME this does not work
  cb_done(h);
elseif isequal(eventdata.Key, 'n') && any(strcmp(eventdata.Modifier, 'shift')) && any(strcmp(eventdata.Modifier, 'control'))
  % FIXME this does not always work correctly
  cb_skip(h);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_show(h, eventdata)
h = parentfig(h);
data = guidata(h);
clf(h);
if data.current<1
  % introduction, construct the GUI elements
  ui_next = uicontrol(h, 'tag', 'ui_next', 'KeyPressFcn', @cb_keyboard, 'style', 'pushbutton', 'string', 'next >', 'callback', @cb_next);
  ui_help = uicontrol(h, 'tag', 'ui_help', 'KeyPressFcn', @cb_keyboard, 'style', 'text', 'max', 10, 'horizontalAlignment', 'left', 'FontName', '', 'backgroundColor', get(h, 'color'));
  % display the introduction text
  title = cat(2, data.file, ' - introduction');
  tmp = {data.script.title};
  help = sprintf('%s\n', tmp{:});
  help = sprintf('This wizard will help you through the following steps:\n\n%s', help);
  set(ui_help, 'string', help);
  set(h, 'Name', title);
elseif data.current>length(data.script)
  % finalization, construct the GUI elements
  ui_prev = uicontrol(h, 'tag', 'ui_prev', 'KeyPressFcn', @cb_keyboard, 'style', 'pushbutton', 'string', '< prev', 'callback', @cb_prev);
  ui_next = uicontrol(h, 'tag', 'ui_next', 'KeyPressFcn', @cb_keyboard, 'style', 'pushbutton', 'string', 'finish', 'callback', @cb_done);
  ui_help = uicontrol(h, 'tag', 'ui_help', 'KeyPressFcn', @cb_keyboard, 'style', 'text', 'max', 10, 'horizontalAlignment', 'left', 'FontName', '', 'backgroundColor', get(h, 'color'));
  % display the finalization text
  title = cat(2, data.file, ' - finish');
  help = sprintf('If you click finish, the variables created by the script will be exported to the MATLAB workspace\n');
  set(ui_help, 'string', help);
  set(h, 'Name', title);
else
  % normal wizard step, construct the GUI elements
  ui_prev = uicontrol(h, 'tag', 'ui_prev', 'KeyPressFcn', @cb_keyboard, 'style', 'pushbutton', 'string', '< prev', 'callback', @cb_prev);
  ui_next = uicontrol(h, 'tag', 'ui_next', 'KeyPressFcn', @cb_keyboard, 'style', 'pushbutton', 'string', 'next >', 'callback', @cb_next);
  ui_help = uicontrol(h, 'tag', 'ui_help', 'KeyPressFcn', @cb_keyboard, 'style', 'text', 'max', 10, 'horizontalAlignment', 'left', 'FontName', '', 'backgroundColor', get(h, 'color'));
  ui_code = uicontrol(h, 'tag', 'ui_code', 'KeyPressFcn', @cb_keyboard, 'style', 'edit', 'max', 10, 'horizontalAlignment', 'left', 'FontName', 'Courier', 'backgroundColor', 'white');
  % display the current wizard content
  help  = data.script(data.current).help;
  code  = data.script(data.current).code;
  title = cat(2, data.file, ' - ', data.script(data.current).title);
  set(ui_help, 'string', help);
  set(ui_code, 'string', code);
  set(h, 'Name', title);
end
cb_resize(h);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_resize(h, eventdata)
h = parentfig(h);
ui_prev = findobj(h, 'tag', 'ui_prev');
ui_next = findobj(h, 'tag', 'ui_next');
ui_help = findobj(h, 'tag', 'ui_help');
ui_code = findobj(h, 'tag', 'ui_code');
siz = get(h, 'position');
x = siz(3);
y = siz(4);
w = (y-40-20)/2;
if ~isempty(ui_prev)
  set(ui_prev, 'position', [x-150 10 60 20]);
end
if ~isempty(ui_next)
  set(ui_next, 'position', [x-080 10 60 20]);
end
if ~isempty(ui_code) && ~isempty(ui_help)
  set(ui_code, 'position', [10 40   x-20 w]);
  set(ui_help, 'position', [10 50+w x-20 w]);
else
  set(ui_help, 'position', [10 40 x-20 y-50]);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_prev(h, eventdata)
h = parentfig(h);
cb_disable(h);
data = guidata(h);
if data.current>0
  data.current = data.current - 1;
end
guidata(h, data);
cb_show(h);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_next(h, eventdata)
h = parentfig(h);
cb_disable(h);
data = guidata(h);
if data.current>0
  title = cat(2, data.file, ' - busy');
  set(h, 'Name', title);
  ui_code = findobj(h, 'tag', 'ui_code');
  code = get(ui_code, 'string');
  try
    for i=1:size(code,1)
      evalin('caller', code(i,:));
    end
  catch
    lasterr;
  end
  data.script(data.current).code = code;
end
data.current = data.current+1;
guidata(h, data);
cb_show(h);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_skip(h, eventdata)
h = parentfig(h);
cb_disable(h);
data = guidata(h);
data.current = data.current+1;
guidata(h, data);
cb_show(h);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_disable(h, eventdata)
h = parentfig(h);
ui_prev = findobj(h, 'tag', 'ui_prev');
ui_next = findobj(h, 'tag', 'ui_next');
ui_help = findobj(h, 'tag', 'ui_help');
ui_code = findobj(h, 'tag', 'ui_code');
set(ui_prev, 'Enable', 'off');
set(ui_next, 'Enable', 'off');
set(ui_help, 'Enable', 'off');
set(ui_code, 'Enable', 'off');
drawnow
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_enable(h, eventdata)
h = parentfig(h);
ui_prev = findobj(h, 'tag', 'ui_prev');
ui_next = findobj(h, 'tag', 'ui_next');
ui_help = findobj(h, 'tag', 'ui_help');
ui_code = findobj(h, 'tag', 'ui_code');
set(ui_prev, 'Enable', 'on');
set(ui_next, 'Enable', 'on');
set(ui_help, 'Enable', 'on');
set(ui_code, 'Enable', 'on');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_edit(h, eventdata)
h = parentfig(h);
data = guidata(h);
filename = fullfile(data.path, [data.file data.ext]);
edit(filename);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_load(h, eventdata)
h = parentfig(h);
cb_disable(h);
[f, p] = uigetfile('*.m', 'Load script from an M-file');
if ~isequal(f,0)
  data = guidata(h);
  filename = fullfile(p, f);
  [p, f, x] = fileparts(filename);
  str = script_parse(filename);
  data.script   = str;
  data.current  = 0;
  data.path = p;
  data.file = f;
  data.ext  = x;
  guidata(h, data);
  cb_show(h);
end
cb_enable(h);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_save(h, eventdata)
cb_disable(h);
data = guidata(h);
filename = fullfile(data.path, [data.file data.ext]);
[f, p] = uiputfile('*.m', 'Save script to an M-file', filename);
if ~isequal(f,0)
  filename = fullfile(p, f);
  [p, f, x] = fileparts(filename);
  fid = fopen(filename, 'wt');
  script = data.script;
  for k=1:length(script)
    fprintf(fid, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    for i=1:size(script(k).help, 1)
      fprintf(fid, '%% %s\n', deblank(script(k).help(i,:)));
    end
    for i=1:size(script(k).code, 1)
      fprintf(fid, '%s\n', deblank(script(k).code(i,:)));
    end
  end
  fclose(fid);
  data.path = p;
  data.file = f;
  data.ext  = x;
  guidata(h, data);
  cb_show(h);
end
cb_enable(h);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_done(h, eventdata)
h = parentfig(h);
cb_disable(h);
assignin('caller', 'wizard_ok', 1);
delete(h);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_cancel(h, eventdata)
h = parentfig(h);
cb_disable(h);
assignin('caller', 'wizard_ok', 0);
delete(h);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function script = script_parse(filename)
if isempty(filename)
  script.title = '';
  script.help  = '';
  script.code  = '';
  return
end
fid = fopen(filename);
str = {};
while ~feof(fid)
  str{end+1} = fgetl(fid);
end
str = str(:);
fclose(fid);
i = 1; % line number
k = 1; % block number
script = [];
while i<=length(str)
  script(k).title = {};
  script(k).help  = {};
  script(k).code = {};
  while i<=length(str) && ~isempty(regexp(str{i}, '^%%%%', 'once'))
    % skip the seperator lines
    i = i+1;
  end
  while i<=length(str) && ~isempty(regexp(str{i}, '^%', 'once'))
    script(k).help{end+1} = str{i}(2:end);
    i = i+1;
  end
  while i<=length(str) && isempty(regexp(str{i}, '^%%%%', 'once'))
    script(k).code{end+1} = str{i};
    i = i+1;
  end
  k = k+1;
end
sel = false(size(script));
for k=1:length(script)
  sel(k) = isempty(script(k).help) && isempty(script(k).code);
  if length(script(k).help)>0
    script(k).title = char(script(k).help{1});
  else
    script(k).title = '<unknown>';
  end
  script(k).help  = char(script(k).help(:));
  script(k).code  = char(script(k).code(:));
end
script(sel) = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = parentfig(h)
while get(h, 'parent')
  h = get(h, 'parent');
end
return
