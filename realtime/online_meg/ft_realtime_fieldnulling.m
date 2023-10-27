function ft_realtime_fieldnulling(cfg)

% FT_REALTIME_FIELDNULLING is a real-time application to drive the nulling
% coils in the magnetically shielded room.
%
% Use as
%   ft_realtime_fieldnulling(cfg)
%
% The configuration should contain
%   cfg.serialport  = string, name of the serial port (default = 'COM2')
%   cfg.fsample     = sampling frequency (default = 1000)

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% set the defaults
cfg.serialport  = ft_getopt(cfg, 'serialport', 'COM2');
cfg.baudrate    = ft_getopt(cfg, 'baudrate', 921600);
cfg.fsample     = ft_getopt(cfg, 'fsample', 400);

% these are not used at the moment, the defaults seem to work fine
cfg.databits    = 8;
cfg.flowcontrol = 'none';
cfg.stopbits    = 1;
cfg.parity      = 'none';

%% set up the serial connection to the fluxgate sensor

% % open the serial device
% fluxgate = serialport(cfg.serialport, cfg.baudrate);
% cleanup = onCleanup(@()serial_cleanup(fluxgate));
%
% %  make sure the persistent variables are not reused from the last call
% clear serial_callback
%
% configureTerminator(fluxgate, 'LF');
% configureCallback(fluxgate, 'terminator', @serial_callback);

%% set up the digital-to-analog converter

% % Create a DataAcquisition object for the specified vendor.
% dac = daq('ni');
%
% % Add channels and set channel properties, if any.
% addoutput(dac,'cDAQ1Mod1','ao0','Voltage');
% addoutput(dac,'cDAQ1Mod1','ao1','Voltage');
% addoutput(dac,'cDAQ1Mod1','ao2','Voltage');
% addoutput(dac,'cDAQ1Mod1','ao3','Voltage');
% addoutput(dac,'cDAQ1Mod1','ao4','Voltage');
% addoutput(dac,'cDAQ1Mod1','ao5','Voltage');
%
% % Output the specified DC amplitude on each channel.
% dcOutput = [1 2 3 4 5 6];
% write(dac, dcOutput);

%%

% open a new figure with the specified settings
fig = open_figure(keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'}));

% store the configuration details and in the application
setappdata(fig, 'cfg', cfg);
setappdata(fig, 'current', [0 0 0 0 0 0 ]);

% add the callbacks
set(fig, 'CloseRequestFcn',     @cb_quit);
set(fig, 'WindowKeypressfcn',   @cb_keyboard);
set(fig, 'WindowButtondownfcn', @cb_click);

cb_creategui(fig);
cb_redraw(fig);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble provenance

if ~ft_nargout
  % don't return anything
  clear cfg
end

end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function serial_callback(fluxgate, varargin)

% this is called on every line/sample that is received
% but the data is not written that fast, so we need a buffer
persistent counter

if isempty(counter)
  counter = 0;
end

% The format is as follows:
% 1353411242214;0.0010416524;-0.0010926605;0.0014391395;0.0020856819<CR><LF>
%
% The separate values have the following meanings:
% Timestamp;Value Channel 1;Value Channel 2;Value Channel 3;Absolute Value<CR><LF>

line = char(readline(fluxgate));
line = strrep(line, ',', '.'); % depending on the language settings there might be a . or a , as decimal separator
dat = str2double(split(line, ';'));

counter = counter + 1;

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function serial_cleanup(fluxgate)
disp('cleanup');
configureCallback(fluxgate, 'off');
disp('stopped')
clear fluxgate
clear readSerialData
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_creategui(h, eventdata, handles)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

p1 = uipanel(fig, 'units', 'normalized', 'position', [0.0 0.5 0.5 0.5]);
p2 = uipanel(fig, 'units', 'normalized', 'position', [0.5 0.5 0.5 0.5]);
p3 = uipanel(fig, 'units', 'normalized', 'position', [0.0 0.0 1.0 0.5]);

uicontrol('parent', p1, 'style', 'text', 'string', 'current', 'units', 'normalized', 'position', [0.0 0.9 0.5 0.1]);
uicontrol('parent', p2, 'style', 'text', 'string', 'field',   'units', 'normalized', 'position', [0.0 0.9 0.5 0.1]);

uicontrol('parent', p1, 'tag', 'c1l', 'style', 'text', 'string', 'x');
uicontrol('parent', p1, 'tag', 'c1e', 'style', 'edit', 'string', '0');
uicontrol('parent', p1, 'tag', 'c1m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c1p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c2e', 'style', 'edit', 'string', '0');
uicontrol('parent', p1, 'tag', 'c2p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c2m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c2x', 'style', 'checkbox');

uicontrol('parent', p1, 'tag', 'c3l', 'style', 'text', 'string', 'y');
uicontrol('parent', p1, 'tag', 'c3e', 'style', 'edit', 'string', '0');
uicontrol('parent', p1, 'tag', 'c3m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c3p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c4e', 'style', 'edit', 'string', '0');
uicontrol('parent', p1, 'tag', 'c4m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c4p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c4x', 'style', 'checkbox');

uicontrol('parent', p1, 'tag', 'c5l', 'style', 'text', 'string', 'z');
uicontrol('parent', p1, 'tag', 'c5e', 'style', 'edit', 'string', '0');
uicontrol('parent', p1, 'tag', 'c5m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c5p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c6e', 'style', 'edit', 'string', '0');
uicontrol('parent', p1, 'tag', 'c6m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c6p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c6x', 'style', 'checkbox');

ft_uilayout(p1, 'style', 'edit', 'callback', @cb_click);
ft_uilayout(p1, 'style', 'checkbox', 'callback', @cb_click);
ft_uilayout(p1, 'style', 'pushbutton', 'callback', @cb_click);

ft_uilayout(p1, 'tag', 'c..', 'units', 'normalized', 'height', 0.1)
ft_uilayout(p1, 'tag', 'c..', 'style', 'edit', 'width', 0.20);
ft_uilayout(p1, 'tag', 'c.p', 'width', 0.05)
ft_uilayout(p1, 'tag', 'c.m', 'width', 0.05)
ft_uilayout(p1, 'tag', 'c.x', 'width', 0.05)

ft_uilayout(p1, 'tag', 'c[12].', 'hpos', 'auto', 'vpos', 0.7);
ft_uilayout(p1, 'tag', 'c[34].', 'hpos', 'auto', 'vpos', 0.5);
ft_uilayout(p1, 'tag', 'c[56].', 'hpos', 'auto', 'vpos', 0.3);

uicontrol('parent', p2, 'tag', 'f1l', 'style', 'text', 'string', 'x');
uicontrol('parent', p2, 'tag', 'f1e', 'style', 'edit', 'string', '0');
uicontrol('parent', p2, 'tag', 'f2l', 'style', 'text', 'string', 'y');
uicontrol('parent', p2, 'tag', 'f2e', 'style', 'edit', 'string', '0');
uicontrol('parent', p2, 'tag', 'f3l', 'style', 'text', 'string', 'x');
uicontrol('parent', p2, 'tag', 'f3e', 'style', 'edit', 'string', '0');
uicontrol('parent', p2, 'tag', 'f4l', 'style', 'text', 'string', 'abs');
uicontrol('parent', p2, 'tag', 'f4e', 'style', 'edit', 'string', '0');

ft_uilayout(p2, 'tag', 'f..', 'units', 'normalized', 'height', 0.1)
ft_uilayout(p2, 'tag', 'f1.', 'hpos', 'auto', 'vpos', 0.7);
ft_uilayout(p2, 'tag', 'f2.', 'hpos', 'auto', 'vpos', 0.5);
ft_uilayout(p2, 'tag', 'f3.', 'hpos', 'auto', 'vpos', 0.3);
ft_uilayout(p2, 'tag', 'f4.', 'hpos', 'auto', 'vpos', 0.1);

uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'calibrate', 'callback', @cb_click);
uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'auto null', 'callback', @cb_click);
uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'coils off', 'callback', @cb_click);
uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'quit','callback', @cb_quit);

ft_uilayout(p3, 'style', 'pushbutton', 'units', 'normalized', 'height', 0.2, 'width', 0.3, 'hpos', 0.35, 'vpos', 'auto');

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

current = getappdata(fig, 'current');
set(findall(fig, 'tag', 'c1e'), 'string', current(1));
set(findall(fig, 'tag', 'c2e'), 'string', current(2));
set(findall(fig, 'tag', 'c3e'), 'string', current(3));
set(findall(fig, 'tag', 'c4e'), 'string', current(4));
set(findall(fig, 'tag', 'c5e'), 'string', current(5));
set(findall(fig, 'tag', 'c6e'), 'string', current(6));

uiresume;
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_click(h, eventdata)
fig = getparent(h);
current = getappdata(fig, 'current');
step = 0.1;
switch get(h, 'tag')
  case 'c1e'
    current(1) = str2double(get(h, 'string'));
  case 'c2e'
    current(2) = str2double(get(h, 'string'));
  case 'c3e'
    current(3) = str2double(get(h, 'string'));
  case 'c4e'
    current(4) = str2double(get(h, 'string'));
  case 'c5e'
    current(5) = str2double(get(h, 'string'));
  case 'c6e'
    current(6) = str2double(get(h, 'string'));
  case 'c1p'
    current(1) = current(1) + step;
  case 'c1m'
    current(1) = current(1) - step;
  case 'c2p'
    current(2) = current(2) + step;
  case 'c2m'
    current(2) = current(2) - step;
  case 'c3p'
    current(3) = current(3) + step;
  case 'c3m'
    current(3) = current(3) - step;
  case 'c4p'
    current(4) = current(4) + step;
  case 'c4m'
    current(4) = current(4) - step;
  case 'c5p'
    current(5) = current(5) + step;
  case 'c5m'
    current(5) = current(5) - step;
  case 'c6p'
    current(6) = current(6) + step;
  case 'c6m'
    current(6) = current(6) - step;
  case 'c2x'
    if get(h, 'Value')
      set(findall(fig, 'tag', 'c2e'), 'Enable', 'off');
      set(findall(fig, 'tag', 'c2m'), 'Enable', 'off');
      set(findall(fig, 'tag', 'c2p'), 'Enable', 'off');
    else
      set(findall(fig, 'tag', 'c2e'), 'Enable', 'on');
      set(findall(fig, 'tag', 'c2m'), 'Enable', 'on');
      set(findall(fig, 'tag', 'c2p'), 'Enable', 'on');
    end
  case 'c4x'
    if get(h, 'Value')
      set(findall(fig, 'tag', 'c4e'), 'Enable', 'off');
      set(findall(fig, 'tag', 'c4m'), 'Enable', 'off');
      set(findall(fig, 'tag', 'c4p'), 'Enable', 'off');
    else
      set(findall(fig, 'tag', 'c4e'), 'Enable', 'on');
      set(findall(fig, 'tag', 'c4m'), 'Enable', 'on');
      set(findall(fig, 'tag', 'c4p'), 'Enable', 'on');
    end
  case 'c6x'
    if get(h, 'Value')
      set(findall(fig, 'tag', 'c6e'), 'Enable', 'off');
      set(findall(fig, 'tag', 'c6m'), 'Enable', 'off');
      set(findall(fig, 'tag', 'c6p'), 'Enable', 'off');
    else
      set(findall(fig, 'tag', 'c6e'), 'Enable', 'on');
      set(findall(fig, 'tag', 'c6m'), 'Enable', 'on');
      set(findall(fig, 'tag', 'c6p'), 'Enable', 'on');
    end
end
if get(findall(fig, 'tag', 'c2x'), 'value')
  current(2) = current(1);
end
if get(findall(fig, 'tag', 'c4x'), 'value')
  current(4) = current(3);
end
if get(findall(fig, 'tag', 'c6x'), 'value')
  current(6) = current(5);
end
setappdata(fig, 'current', current);
cb_redraw(fig);
end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)
fig = getparent(h);

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(fig, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

if isempty(key)
  % this happens if you press the apple key
  key = '';
end

switch key
  case {'' 'shift+shift' 'alt-alt' 'control+control' 'command-0'}
    % do nothing

  case 'q'
    cb_quit(h);

  otherwise
    % do nothing

end % switch key
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_quit(h, eventdata)
fig = getparent(h);
delete(fig)
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end
end % function
