function ft_realtime_fieldnulling(cfg)

% FT_REALTIME_FIELDNULLING is a real-time application to drive the nulling
% coilss in the magnetically shielded room.
%
% Use as
%   ft_realtime_fieldnulling(cfg)
%
% The configuration should contain
%   cfg.serialport  = string, name of the serial port (default = 'COM2')
%   cfg.fsample     = sampling frequency (default = 1000)
%
% When clicking the + and - buttons, you can use the shift and ctrl
% modifier keys to make small and even smaller steps.
%
% See also FT_REALTIME_SENSYSPROXY

% Copyright (C) 2023, Robert Oostenveld
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
cfg.fsample     = ft_getopt(cfg, 'fsample', 40);
cfg.databits    = ft_getopt(cfg, 'databits', 8);
cfg.flowcontrol = ft_getopt(cfg, 'flowcontrol', 'none');
cfg.stopbits    = ft_getopt(cfg, 'stopbits', 1);
cfg.parity      = ft_getopt(cfg, 'parity', 'none');
cfg.position    = ft_getopt(cfg, 'position'); % default is handled below

if isempty(cfg.position)
    cfg.position = get(groot, 'defaultFigurePosition');
    cfg.position(4) = 240;
end

%%

% open a new figure with the specified settings
tmpcfg = keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'});
% overrule some of the settings
tmpcfg.figure = 'ui';
tmpcfg.figurename = 'ft_realtime_fieldnulling';
fig = open_figure(tmpcfg);
drawnow

% store the configuration details and in the application
setappdata(fig, 'cfg', cfg);
setappdata(fig, 'field', [0 0 0 0]); % x, y, z, abs
setappdata(fig, 'offset', [0 0 0 0 0 0]);
setappdata(fig, 'running', false);

% add the callbacks
set(fig, 'CloseRequestFcn',     @cb_quit);
set(fig, 'WindowKeypressfcn',   @cb_keyboard);
set(fig, 'WindowButtondownfcn', @cb_click);

%% set up the serial connection to the fluxgate sensor

measure_callback(fig); % call it once to pass the figure handle

fluxgate = serialport(cfg.serialport, cfg.baudrate);
cleanup = onCleanup(@()measure_cleanup(fluxgate));
configureTerminator(fluxgate, 'LF');
configureCallback(fluxgate, 'terminator', @measure_callback);

%% set up the digital-to-analog converter

coils = daq('ni');
% add channels and set channel properties, if any.
addoutput(coils, 'cDAQ1Mod1', 'ao0', 'Voltage');
addoutput(coils, 'cDAQ1Mod1', 'ao1', 'Voltage');
addoutput(coils, 'cDAQ1Mod1', 'ao2', 'Voltage');
addoutput(coils, 'cDAQ1Mod1', 'ao3', 'Voltage');
addoutput(coils, 'cDAQ1Mod1', 'ao4', 'Voltage');
addoutput(coils, 'cDAQ1Mod1', 'ao5', 'Voltage');

%% activate the graphical user interface

setappdata(fig, 'fluxgate', fluxgate);
setappdata(fig, 'coils', coils);
setappdata(fig, 'running', true);

cb_creategui(fig);
cb_redraw(fig);

while ishandle(fig)
    pause(0.1);
end

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
function control_output_dc(coils, offset)
% output the specified DC amplitude on each channel.
nchan = 6;
assert(numel(offset)==nchan);
write(coils, offset);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function control_output_calibration(h, eventdata)
% output a sine wave with a different frequency on each of the coils
fig = getparent(h);
coils = getappdata(fig, 'coils');

offset = [0 0 0 0 0 0];
frequency = [5 6 7 8 9 10];
amplitude = 1;
duration = 5;
fsample = 1000;
nchan = 6;

assert(numel(offset)==nchan);
assert(numel(frequency)==nchan);
assert(numel(amplitude)==1);
assert(numel(duration)==1);

coils.Rate = fsample;

signal = zeros(nchan,fsample);
time = (0:(fsample-1))/fsample;
for i=1:nchan
    signal(i,:) = amplitude * sin(frequency(i)*2*pi*time) + offset(i);
end
preload(coils, signal');

disp('start calibration');
t = timer('StartDelay', duration, 'ExecutionMode', 'singleshot');
t.StartFcn = @(varargin) start(coils, 'repeatoutput');
t.TimerFcn = @(varargin) stop(coils);
t.StopFcn  = @(varargin) disp('stop calibration');
start(t)

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measure_callback(fluxgate, varargin)

persistent fig counter

if ishandle(fluxgate)
    % it is called once so that the callback knows about the main figure
    fig = fluxgate;
    return
end

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

% add the fluxgate field strength to the figure
setappdata(fig, 'field', dat(2:5));
cb_redraw(fig);

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measure_cleanup(fluxgate)
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

p1 = uipanel(fig, 'units', 'normalized', 'position', [0.0 0.3 0.5 0.7]);
p2 = uipanel(fig, 'units', 'normalized', 'position', [0.5 0.3 0.5 0.7]);
p3 = uipanel(fig, 'units', 'normalized', 'position', [0.0 0.0 1.0 0.3]);

uicontrol('parent', p1, 'style', 'text', 'string', 'control', 'units', 'normalized', 'position', [0.0 0.9 0.5 0.1]);
uicontrol('parent', p2, 'style', 'text', 'string', 'measure',   'units', 'normalized', 'position', [0.0 0.9 0.5 0.1]);

uicontrol('parent', p1, 'tag', 'c1l', 'style', 'text', 'string', 'x');
uicontrol('parent', p1, 'tag', 'c1e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c1m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c1p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c2e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c2p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c2m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c2x', 'style', 'checkbox');

uicontrol('parent', p1, 'tag', 'c3l', 'style', 'text', 'string', 'y');
uicontrol('parent', p1, 'tag', 'c3e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c3m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c3p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c4e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c4m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c4p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c4x', 'style', 'checkbox');

uicontrol('parent', p1, 'tag', 'c5l', 'style', 'text', 'string', 'z');
uicontrol('parent', p1, 'tag', 'c5e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c5m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c5p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c6e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c6m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c6p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c6x', 'style', 'checkbox');

ft_uilayout(p1, 'style', 'edit', 'callback', @cb_click, 'backgroundcolor', [1 1 1]);
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
uicontrol('parent', p2, 'tag', 'f1e', 'style', 'edit');
uicontrol('parent', p2, 'tag', 'f2l', 'style', 'text', 'string', 'y');
uicontrol('parent', p2, 'tag', 'f2e', 'style', 'edit');
uicontrol('parent', p2, 'tag', 'f3l', 'style', 'text', 'string', 'x');
uicontrol('parent', p2, 'tag', 'f3e', 'style', 'edit');
uicontrol('parent', p2, 'tag', 'f4l', 'style', 'text', 'string', 'abs');
uicontrol('parent', p2, 'tag', 'f4e', 'style', 'edit');

ft_uilayout(p2, 'tag', 'f..', 'units', 'normalized', 'height', 0.1)
ft_uilayout(p2, 'tag', 'f1.', 'hpos', 'auto', 'vpos', 0.7);
ft_uilayout(p2, 'tag', 'f2.', 'hpos', 'auto', 'vpos', 0.5);
ft_uilayout(p2, 'tag', 'f3.', 'hpos', 'auto', 'vpos', 0.3);
ft_uilayout(p2, 'tag', 'f4.', 'hpos', 'auto', 'vpos', 0.1);

uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'calibrate', 'callback', @control_output_calibration);
uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'auto null');
uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'quit',      'callback', @cb_quit);

ft_uilayout(p3, 'style', 'pushbutton', 'units', 'normalized', 'width', 0.35, 'hpos', 'auto', 'vpos', 0.5);

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)
fig = getparent(h);
offset = getappdata(fig, 'offset');
field = getappdata(fig, 'field');
coils = getappdata(fig, 'coils');

if get(findall(fig, 'tag', 'c2x'), 'value')
    offset(2) = offset(1);
end
if get(findall(fig, 'tag', 'c4x'), 'value')
    offset(4) = offset(3);
end
if get(findall(fig, 'tag', 'c6x'), 'value')
    offset(6) = offset(5);
end

if getappdata(fig, 'running')
    control_output_dc(coils, offset);

    % update the current driver offset
    set(findall(fig, 'tag', 'c1e'), 'string', offset(1));
    set(findall(fig, 'tag', 'c2e'), 'string', offset(2));
    set(findall(fig, 'tag', 'c3e'), 'string', offset(3));
    set(findall(fig, 'tag', 'c4e'), 'string', offset(4));
    set(findall(fig, 'tag', 'c5e'), 'string', offset(5));
    set(findall(fig, 'tag', 'c6e'), 'string', offset(6));

    % update the fluxgate field
    set(findall(fig, 'tag', 'f1e'), 'string', sprintf('%.03f uT', 1e6*field(1)));
    set(findall(fig, 'tag', 'f2e'), 'string', sprintf('%.03f uT', 1e6*field(2)));
    set(findall(fig, 'tag', 'f3e'), 'string', sprintf('%.03f uT', 1e6*field(3)));
    set(findall(fig, 'tag', 'f4e'), 'string', sprintf('%.03f uT', 1e6*field(4)));
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_click(h, eventdata)
fig = getparent(h);
offset = getappdata(fig, 'offset');
switch get(fig, 'SelectionType')
    case 'alt'       % ctrl-click
        step = 0.01;
    case 'extend'    % shift-click
        step = 0.1;
    otherwise
        step = 1;
end
switch get(h, 'tag')
    case 'c1e'
        offset(1) = str2double(get(h, 'string'));
    case 'c2e'
        offset(2) = str2double(get(h, 'string'));
    case 'c3e'
        offset(3) = str2double(get(h, 'string'));
    case 'c4e'
        offset(4) = str2double(get(h, 'string'));
    case 'c5e'
        offset(5) = str2double(get(h, 'string'));
    case 'c6e'
        offset(6) = str2double(get(h, 'string'));
    case 'c1p'
        offset(1) = offset(1) + step;
    case 'c1m'
        offset(1) = offset(1) - step;
    case 'c2p'
        offset(2) = offset(2) + step;
    case 'c2m'
        offset(2) = offset(2) - step;
    case 'c3p'
        offset(3) = offset(3) + step;
    case 'c3m'
        offset(3) = offset(3) - step;
    case 'c4p'
        offset(4) = offset(4) + step;
    case 'c4m'
        offset(4) = offset(4) - step;
    case 'c5p'
        offset(5) = offset(5) + step;
    case 'c5m'
        offset(5) = offset(5) - step;
    case 'c6p'
        offset(6) = offset(6) + step;
    case 'c6m'
        offset(6) = offset(6) - step;
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

% set near-zero values to zero
offset(abs(offset)<10*eps) = 0;

% clip to the allowed range
offset(offset>+10) = +10;
offset(offset<-10) = -10;

setappdata(fig, 'offset', offset);
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
setappdata(fig, 'running', false);
delete(fig);
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
