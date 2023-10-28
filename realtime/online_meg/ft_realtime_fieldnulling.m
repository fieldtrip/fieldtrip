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
cfg.serialport    = ft_getopt(cfg, 'serialport', 'COM2');
cfg.baudrate      = ft_getopt(cfg, 'baudrate', 921600);
cfg.fsample       = ft_getopt(cfg, 'fsample', 40);
cfg.databits      = ft_getopt(cfg, 'databits', 8);
cfg.flowcontrol   = ft_getopt(cfg, 'flowcontrol', 'none');
cfg.stopbits      = ft_getopt(cfg, 'stopbits', 1);
cfg.parity        = ft_getopt(cfg, 'parity', 'none');
cfg.position      = ft_getopt(cfg, 'position'); % default is handled below
cfg.enableinput   = ft_getopt(cfg, 'enableinput', 'yes');
cfg.enableoutput  = ft_getopt(cfg, 'enableoutput', 'yes');
cfg.offset        = ft_getopt(cfg, 'offset', [0 0 0 0 0 0]); % initial offset on the coils

% these determine the calibration signal
cfg.calibration           = ft_getopt(cfg, '', []);
cfg.calibration.duration  = ft_getopt(cfg, 'duration', 5); % in seconds
cfg.calibration.nchan     = ft_getopt(cfg, 'nchan', 6);
cfg.calibration.amplitude = ft_getopt(cfg, 'amplitude', 1);
cfg.calibration.frequency = ft_getopt(cfg, 'frequency', [5 6 7 8 9 10]); % in Hz
cfg.calibration.fsample   = ft_getopt(cfg, 'fsample', 1000); % in Hz

if isempty(cfg.position)
  cfg.position = get(groot, 'defaultFigurePosition');
  cfg.position(4) = 240;
end

% open a new figure with the specified settings
tmpcfg = keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'});
% overrule some of the settings
tmpcfg.figure = 'ui';
tmpcfg.figurename = 'ft_realtime_fieldnulling';
fig = open_figure(tmpcfg);
drawnow

% store the configuration details in the application
setappdata(fig, 'cfg', cfg);
setappdata(fig, 'field', [nan nan nan nan]);  % set the initial measured field: x, y, z, abs
setappdata(fig, 'offset', cfg.offset);        % set the initial offset
setappdata(fig, 'running', false);
setappdata(fig, 'calib', []);
setappdata(fig, 'calibration', 'off');

% add the callbacks
set(fig, 'CloseRequestFcn',     @cb_quit);
set(fig, 'WindowKeypressfcn',   @cb_keyboard);
set(fig, 'WindowButtondownfcn', @cb_click);

%% set up the serial connection to the fluxgate sensor

measure_callback(fig); % call it once to pass the figure handle

if istrue(cfg.enableinput)
  fluxgate = serialport(cfg.serialport, cfg.baudrate);
  cleanup = onCleanup(@()measure_cleanup(fluxgate));
  configureTerminator(fluxgate, 'LF');
  configureCallback(fluxgate, 'terminator', @measure_callback);
else
  fluxgate = [];
end

setappdata(fig, 'fluxgate', fluxgate);

%% set up the digital-to-analog converter

if istrue(cfg.enableoutput)
  coils = daq('ni');
  % add channels and set channel properties, if any.
  addoutput(coils, 'cDAQ1Mod1', 'ao0', 'Voltage');
  addoutput(coils, 'cDAQ1Mod1', 'ao1', 'Voltage');
  addoutput(coils, 'cDAQ1Mod1', 'ao2', 'Voltage');
  addoutput(coils, 'cDAQ1Mod1', 'ao3', 'Voltage');
  addoutput(coils, 'cDAQ1Mod1', 'ao4', 'Voltage');
  addoutput(coils, 'cDAQ1Mod1', 'ao5', 'Voltage');
else
  coils = [];
end

setappdata(fig, 'coils', coils);

%% activate the graphical user interface

setappdata(fig, 'running', true);
cb_create_gui(fig);
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
function control_offset(fig)
coils = getappdata(fig, 'coils');
offset = getappdata(fig, 'offset');

if ~isempty(coils)
  % output the specified DC offset on each of the coils
  write(coils, offset);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function control_auto_null(h, eventdata)
fig = getparent(h);
offset = getappdata(fig, 'offset');
field = getappdata(fig, 'field');
calib = getappdata(fig, 'calib');

% only keep the x, y, and z component
residual = field(1:3);
residual = residual(:); % it should be a column
offset = offset(:); % it should be a column

if isempty(calib)
  warning('cannot auto-null, calibration has not yet been performed')

else
  disp('auto-null')

  % The idea is to compute a correction to the current offset that will bring the measured field to zero
  %   field = residual - calib * offset  % this is the environmental field
  %   residual = field + calib * offset  % this is the residual field that we measure
  %
  % The sum of the enviromental field, the offset, and the correction should be zero
  %   zero = field + calib * offset + calib * correction
  %        = residual               + calib * correction
  % Hence
  %   residual = - calib * correction

  correction = - calib \ residual;
  offset = offset + correction;

  % show and apply the updated offset values
  setappdata(fig, 'offset', offset);
  cb_redraw(fig);
  control_offset(fig);

end
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibration_init(h, eventdata)
fig = getparent(h);
cb_disable_gui(fig);
setappdata(fig', 'calibration', 'init');
disp('initiated calibration')
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibration_start_signal(fig)
% output a sine wave with a different frequency on each of the coils
cfg     = getappdata(fig, 'cfg');
coils   = getappdata(fig, 'coils');
offset  = getappdata(fig, 'offset');
nchan   = cfg.calibration.nchan;
fsample = cfg.calibration.fsample;

assert(numel(offset)==nchan);
assert(numel(cfg.calibration.frequency)==nchan);
assert(numel(cfg.calibration.amplitude)==1);

coils.Rate = fsample;

% create one second of data, it will loop until finished
signal = zeros(nchan,fsample);
time = (0:(fsample-1))/fsample;
for i=1:nchan
  signal(i,:) = cfg.calibration.amplitude * sin(cfg.calibration.frequency(i)*2*pi*time) + offset(i);
end

if ~isempty(coils)
  preload(coils, signal');
  start(coils);
  disp('started calibration signal');
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibration_stop_signal(fig)
coils = getappdata(fig, 'coils');
if ~isempty(coils)
  stop(coils);
  disp('stopped calibration signal');
end
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibration_compute(fig, dat)
cfg     = getappdata(fig, 'cfg');
nchan   = cfg.calibration.nchan;
fsample = cfg.calibration.fsample;
nsample = size(dat,2);

time = (0:(nsample-1))/fsample;
input = zeros(nchan+1,nsample); % additional channel for the constant offset
for i=1:nchan
  input(i,:) = sin(cfg.calibration.frequency(i)*2*pi*time);
end
input(end,:) = 1;

% compute the linear mix from the input signals on the coils towards the output measurement on the fluxgate
% output = calib * input + noise
calib = dat / input;

disp(calib)

noise = dat - calib*input;
gof = 1 - norm(noise, 'fro')/norm(dat, 'fro');

fprintf('calibration computed\n');
fprintf('goodness of fit = %f %%\n', gof*100);

setappdata(fig, 'calib', calib);
cb_enable_gui(fig);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measure_callback(fluxgate, varargin)
persistent fig counter buffer

if ishandle(fluxgate)
  % it is called like this once, so that it knows about the main figure
  fig = fluxgate;
  return
end

% The serial communication format is as follows:
% 1353411242214;0.0010416524;-0.0010926605;0.0014391395;0.0020856819<CR><LF>
%
% The separate values have the following meanings:
% Timestamp;Value Channel 1;Value Channel 2;Value Channel 3;Absolute Value<CR><LF>

line = char(readline(fluxgate));
line = strrep(line, ',', '.'); % depending on the language settings there might be a . or a , as decimal separator
dat = str2double(split(line, ';'));

% the calibration is implemented as a finite-state machine 
% it switches from 'init' -> 'on' -> 'off'
switch getappdata(fig, 'calibration')
  case 'init'
    setappdata(fig, 'calibration', 'on');
    calibration_start_signal(fig);
    cfg = getappdata(fig, 'cfg');
    nsample = round(cfg.calibration.duration * cfg.fsample);
    nchan = numel(dat);
    buffer = nan(nchan, nsample);
    counter = 0;
  case 'on'
    if counter<nsamples
      % add the current sample to the buffer
      counter = counter+1;
      buffer(:,counter) = dat;
    else
      % compute the calibration values and cleanup
      setappdata(fig, 'calibration', 'off');
      calibration_stop_signal(fig);
      % the first channel is the timestamp, the last one is the absolute value
      dat = buffer([2 3 4],:);
      calibration_compute(fig, dat);
      buffer = [];
      counter = 0;
    end
  case 'off'
    % nothing to do
end

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
function cb_create_gui(h, eventdata, handles)
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

ft_uilayout(p1, 'tag', 'c..', 'units', 'normalized', 'height', 0.15)
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

uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'calibrate', 'callback', @calibration_init);
uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'auto null', 'callback', @control_auto_null);
uicontrol('parent', p3, 'style', 'pushbutton', 'string', 'quit',      'callback', @cb_quit);

ft_uilayout(p3, 'style', 'pushbutton', 'units', 'normalized', 'width', 0.35, 'height', 0.4, 'hpos', 'auto', 'vpos', 0.5);

cb_enable_gui(fig);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_enable_gui(h, eventdata)
disp('enable gui');
fig = getparent(h);
set(findall(fig, 'type', 'uicontrol', 'style', 'edit'),       'Enable', 'on');
set(findall(fig, 'type', 'uicontrol', 'style', 'pushbutton'), 'Enable', 'on');
set(findall(fig, 'type', 'uicontrol', 'style', 'checkbox'),   'Enable', 'on');

h = findall(fig, 'tag', 'c2x');
if get(h, 'Value')
  set(findall(fig, 'tag', 'c2e'), 'Enable', 'off');
  set(findall(fig, 'tag', 'c2m'), 'Enable', 'off');
  set(findall(fig, 'tag', 'c2p'), 'Enable', 'off');
else
  set(findall(fig, 'tag', 'c2e'), 'Enable', 'on');
  set(findall(fig, 'tag', 'c2m'), 'Enable', 'on');
  set(findall(fig, 'tag', 'c2p'), 'Enable', 'on');
end

h = findall(fig, 'tag', 'c4x');
if get(h, 'Value')
  set(findall(fig, 'tag', 'c4e'), 'Enable', 'off');
  set(findall(fig, 'tag', 'c4m'), 'Enable', 'off');
  set(findall(fig, 'tag', 'c4p'), 'Enable', 'off');
else
  set(findall(fig, 'tag', 'c4e'), 'Enable', 'on');
  set(findall(fig, 'tag', 'c4m'), 'Enable', 'on');
  set(findall(fig, 'tag', 'c4p'), 'Enable', 'on');
end

h = findall(fig, 'tag', 'c6x');
if get(h, 'Value')
  set(findall(fig, 'tag', 'c6e'), 'Enable', 'off');
  set(findall(fig, 'tag', 'c6m'), 'Enable', 'off');
  set(findall(fig, 'tag', 'c6p'), 'Enable', 'off');
else
  set(findall(fig, 'tag', 'c6e'), 'Enable', 'on');
  set(findall(fig, 'tag', 'c6m'), 'Enable', 'on');
  set(findall(fig, 'tag', 'c6p'), 'Enable', 'on');
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_disable_gui(h, eventdata)
disp('disable gui');
fig = getparent(h);
set(findall(fig, 'type', 'uicontrol', 'style', 'edit'),       'Enable', 'off');
set(findall(fig, 'type', 'uicontrol', 'style', 'pushbutton'), 'Enable', 'off');
set(findall(fig, 'type', 'uicontrol', 'style', 'checkbox'),   'Enable', 'off');
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)
fig     = getparent(h);
offset  = getappdata(fig, 'offset');
field   = getappdata(fig, 'field');

% deal with the two coils being locked to each other or not
if get(findall(fig, 'tag', 'c2x'), 'value')
  offset(2) = offset(1);
end
if get(findall(fig, 'tag', 'c4x'), 'value')
  offset(4) = offset(3);
end
if get(findall(fig, 'tag', 'c6x'), 'value')
  offset(6) = offset(5);
end
setappdata(fig, 'offset', offset);

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

if getappdata(fig, 'running')
  control_offset(fig);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_click(h, eventdata)
fig     = getparent(h);
offset  = getappdata(fig, 'offset');

switch get(fig, 'SelectionType')
  case 'alt'       % ctrl-click, this does not work on macOS
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
    cb_enable_gui(fig); % disable/enable the linked coils
  case 'c4x'
    cb_enable_gui(fig); % disable/enable the linked coils
  case 'c6x'
    cb_enable_gui(fig); % disable/enable the linked coils
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
