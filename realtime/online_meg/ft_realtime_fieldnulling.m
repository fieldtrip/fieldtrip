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
cfg.serialport        = ft_getopt(cfg, 'serialport', 'COM2');
cfg.baudrate          = ft_getopt(cfg, 'baudrate', 921600);
cfg.fsample           = ft_getopt(cfg, 'fsample', 400);

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
% configureTerminator(fluxgate, "LF");
% configureCallback(fluxgate, "terminator", @serial_callback);

%% set up the digital-to-analog converter

% % Create a DataAcquisition object for the specified vendor.
% dac = daq("ni");
% 
% % Add channels and set channel properties, if any.
% addoutput(dac,"cDAQ1Mod1","ao0","Voltage");
% addoutput(dac,"cDAQ1Mod1","ao1","Voltage");
% addoutput(dac,"cDAQ1Mod1","ao2","Voltage");
% addoutput(dac,"cDAQ1Mod1","ao3","Voltage");
% addoutput(dac,"cDAQ1Mod1","ao4","Voltage");
% addoutput(dac,"cDAQ1Mod1","ao5","Voltage");
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
set(fig, 'WindowButtondownfcn', @cb_buttonpress);

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
configureCallback(fluxgate, "off");
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

uicontrol('tag', 'x1l', 'parent', fig, 'style', 'text', 'string', 'x');
uicontrol('tag', 'x1e', 'parent', fig, 'style', 'edit', 'string', '0');
uicontrol('tag', 'x1p', 'parent', fig, 'style', 'pushbutton', 'string', '+', 'callback', @cb_button);
uicontrol('tag', 'x1m', 'parent', fig, 'style', 'pushbutton', 'string', '-', 'callback', @cb_button);
uicontrol('tag', 'x2e', 'parent', fig, 'style', 'edit', 'string', '0');
uicontrol('tag', 'x2p', 'parent', fig, 'style', 'pushbutton', 'string', '+', 'callback', @cb_button);
uicontrol('tag', 'x2m', 'parent', fig, 'style', 'pushbutton', 'string', '-', 'callback', @cb_button);

uicontrol('tag', 'y1l', 'parent', fig, 'style', 'text', 'string', 'y');
uicontrol('tag', 'y1e', 'parent', fig, 'style', 'edit', 'string', '0');
uicontrol('tag', 'y1p', 'parent', fig, 'style', 'pushbutton', 'string', '+', 'callback', @cb_button);
uicontrol('tag', 'y1m', 'parent', fig, 'style', 'pushbutton', 'string', '-', 'callback', @cb_button);

uicontrol('tag', 'y2e', 'parent', fig, 'style', 'edit', 'string', '0');
uicontrol('tag', 'y2p', 'parent', fig, 'style', 'pushbutton', 'string', '+', 'callback', @cb_button);
uicontrol('tag', 'y2m', 'parent', fig, 'style', 'pushbutton', 'string', '-', 'callback', @cb_button);

uicontrol('tag', 'z1l', 'parent', fig, 'style', 'text', 'string', 'z');
uicontrol('tag', 'z1e', 'parent', fig, 'style', 'edit', 'string', '0');
uicontrol('tag', 'z1p', 'parent', fig, 'style', 'pushbutton', 'string', '+', 'callback', @cb_button);
uicontrol('tag', 'z1m', 'parent', fig, 'style', 'pushbutton', 'string', '-', 'callback', @cb_button);

uicontrol('tag', 'z2e', 'parent', fig, 'style', 'edit', 'string', '0');
uicontrol('tag', 'z2p', 'parent', fig, 'style', 'pushbutton', 'string', '+', 'callback', @cb_button);
uicontrol('tag', 'z2m', 'parent', fig, 'style', 'pushbutton', 'string', '-', 'callback', @cb_button);

ft_uilayout(fig, 'tag', '..e', 'units', 'normalized', 'width', 0.20, 'backgroundcolor', [1 1 1])
ft_uilayout(fig, 'tag', '..p', 'units', 'normalized', 'width', 0.05)
ft_uilayout(fig, 'tag', '..m', 'units', 'normalized', 'width', 0.05)

ft_uilayout(fig, 'tag', 'z..', 'units', 'normalized', 'height', 0.05, 'hpos', 'auto', 'vpos', 0.25, 'backgroundcolor', [0.8 0.8 0.8]);
ft_uilayout(fig, 'tag', 'y..', 'units', 'normalized', 'height', 0.05, 'hpos', 'auto', 'vpos', 0.15, 'backgroundcolor', [0.8 0.8 0.8]);
ft_uilayout(fig, 'tag', 'x..', 'units', 'normalized', 'height', 0.05, 'hpos', 'auto', 'vpos', 0.05, 'backgroundcolor', [0.8 0.8 0.8]);

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)
fig = getparent(h);
cfg = getappdata(fig, 'cfg');

% disable all GUI elements
set(findall(fig, 'type', 'UIControl'), 'Enable', 'off')

% FIXME ...

% re-enable all GUI elements
set(findall(fig, 'type', 'UIControl'), 'Enable', 'on')
uiresume;
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_button(h, eventdata)
fig = getparent(h);
current = getappdata(fig, 'current');
step = 0.1;
switch get(h, 'tag')
    case 'x1p'
        current(1) = current(1) + step;
    case 'x1m'
        current(1) = current(1) - step;
    case 'x2p'
        current(2) = current(2) + step;
    case 'x2m'
        current(2) = current(2) - step;
    case 'y1p'
        current(3) = current(3) + step;
    case 'y1m'
        current(3) = current(3) - step;
    case 'y2p'
        current(4) = current(4) + step;
    case 'y2m'
        current(4) = current(4) - step;
    case 'z1p'
        current(5) = current(5) + step;
    case 'z1m'
        current(5) = current(5) - step;
    case 'z2p'
        current(6) = current(6) + step;
    case 'z2m'
        current(6) = current(6) - step;
end
setappdata(fig, 'current', current);
current
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
