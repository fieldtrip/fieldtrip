function ft_realtime_fieldnulling(cfg)

% FT_REALTIME_FIELDNULLING is a real-time application to drive the nulling
% coils in the magnetically shielded room.
%
% Use as
%   ft_realtime_fieldnulling(cfg)
%
% The configuration should contain
%   cfg.serialport   = string, name of the serial port (default = 'COM2')
%   cfg.baudrate     = number (default = 921600)
%   cfg.fsample      = number, sampling rate of the fluxgate AFTER averaging (no default provided)
%   cfg.enableinput  = string, 'yes' or 'no' to enable the fluxgate input (default = 'yes')
%   cfg.enableoutput = string, 'yes' or 'no' to enable the NI-DAQ output (default = 'yes')
%   cfg.range        = [min max], values outside this range will be clipped (default [-10 10])
%   cfg.voltage      = initial voltage on the coils (default = [0 0 0 0 0 0])
%   cfg.polarity     = vector with +1 or -1 values to indicate polarity of coils (default = [1 1 1 1 1 1])
%
% Furthermore, these options determine the calibration signal
%   cfg.calibration.duration  = number in seconds (default = 5)
%   cfg.calibration.pad       = number in seconds (default = 0)
%   cfg.calibration.amplitude = number in Volt (default = 1)
%   cfg.calibration.frequency = number in Hz for each of the coils (default = [1 2 3 4 5 6])
%
% When clicking the + and - buttons, you can use the shift and ctrl
% modifier keys to make small and even smaller steps. The step size depends
% on the range.
%
% See also FT_REALTIME_SENSYSPROXY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some ideas to implement:
%  - dithering of output DAC values to achive
%  - smooth padded transitions, to prevent overshoots
%  - logging of input and output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  % do not continue function execution if something is wrong
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', 'fsample');
cfg = ft_checkconfig(cfg, 'renamed', {'clip', 'range'});

% set the defaults
cfg.serialport    = ft_getopt(cfg, 'serialport', 'COM2');
cfg.baudrate      = ft_getopt(cfg, 'baudrate', 921600);
cfg.fsample       = ft_getopt(cfg, 'fsample'); % this must be specified
cfg.databits      = ft_getopt(cfg, 'databits', 8);
cfg.flowcontrol   = ft_getopt(cfg, 'flowcontrol', 'none');
cfg.stopbits      = ft_getopt(cfg, 'stopbits', 1);
cfg.parity        = ft_getopt(cfg, 'parity', 'none');
cfg.position      = ft_getopt(cfg, 'position'); % default is handled below
cfg.enableinput   = ft_getopt(cfg, 'enableinput', 'yes');
cfg.enableoutput  = ft_getopt(cfg, 'enableoutput', 'yes');
cfg.voltage       = ft_getopt(cfg, 'voltage', [0 0 0 0 0 0]); % initial voltage on each of the coils
cfg.polarity      = ft_getopt(cfg, 'polarity', [1 1 1 1 1 1]); % +1 or -1 for each of the coils coils
cfg.calib         = ft_getopt(cfg, 'calib'); % user-specified calibration values
cfg.range         = ft_getopt(cfg, 'range', [-10 10]); % values outside this range will be clipped

% these determine the calibration signal
cfg.calibration           = ft_getopt(cfg, 'calibration', []);
cfg.calibration.duration  = ft_getopt(cfg.calibration, 'duration', 5); % in seconds
cfg.calibration.pad       = ft_getopt(cfg.calibration, 'pad', 0); % in seconds
cfg.calibration.amplitude = ft_getopt(cfg.calibration, 'amplitude', 1);
cfg.calibration.frequency = ft_getopt(cfg.calibration, 'frequency', [1 2 3 4 5 6]); % in Hz
cfg.calibration.fsample   = ft_getopt(cfg.calibration, 'fsample', 1000); % in Hz

if isempty(cfg.position)
  cfg.position = get(groot, 'defaultFigurePosition');
  cfg.position(3) = 800;
  cfg.position(4) = 240;
end

% clear the subfunctions with persistent values
% FIXME this seems not to work as expected
clear control_voltage
clear measure_sample

% open a new figure with the specified settings
tmpcfg = keepfields(cfg, {'figure', 'position', 'visible', 'renderer', 'figurename', 'title'});
% overrule some of the settings
tmpcfg.figure = 'ui';
tmpcfg.figurename = 'ft_realtime_fieldnulling';
fig = open_figure(tmpcfg);
drawnow

% store the configuration details in the application
setappdata(fig, 'cfg', cfg);
setappdata(fig, 'field', [0 0 0 0]);          % set the initial value for the field: x, y, z, abs
setappdata(fig, 'voltage', cfg.voltage);      % set the initial voltage
setappdata(fig, 'calib', cfg.calib);
setappdata(fig, 'closedloop', 'off');         % can be 'off', 'on'
setappdata(fig, 'calibration', 'off');        % can be 'init', 'on', 'off'

% add the callbacks
set(fig, 'CloseRequestFcn',     @cb_quit);
set(fig, 'WindowKeypressfcn',   @cb_keyboard);
set(fig, 'WindowButtondownfcn', @cb_click);

%% set up the serial connection to the fluxgate sensor

measure_sample(fig); % call it once to pass the figure handle

if istrue(cfg.enableinput)
  ft_info('initializing serial port for fluxgate');
  fluxgate = serialport(cfg.serialport, cfg.baudrate);
  cleanup = onCleanup(@()measure_cleanup(fluxgate));
  configureTerminator(fluxgate, 'LF');
  configureCallback(fluxgate, 'terminator', @measure_sample);
else
  fluxgate = [];
end

setappdata(fig, 'fluxgate', fluxgate);

%% set up the digital-to-analog converter

if istrue(cfg.enableoutput)
  ft_info('initializing NI-DAQ');
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

ft_info('creating gui');
setappdata(fig, 'closedloop', 'off');
cb_create_gui(fig);

t = timer;
t.Period = 1;
t.ExecutionMode = 'fixedRate';
t.TimerFcn = {@cb_redraw, fig};
t.start();

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
function control_voltage(fig)
persistent previous_voltage
cfg     = getappdata(fig, 'cfg');
coils   = getappdata(fig, 'coils');
voltage = getappdata(fig, 'voltage');

% clip the voltage to the specified range
if any(voltage<cfg.range(1) | voltage>cfg.range(2))
  voltage(voltage>cfg.range(2)) = cfg.range(2);
  voltage(voltage<cfg.range(1)) = cfg.range(1);
  setappdata(fig, 'voltage', voltage);
end

% it should be a row vector
voltage = voltage(:)';

if ~isempty(coils) && ~coils.Running && ~isequal(previous_voltage, voltage)
  % output the specified DC voltage on the coils
  write(coils, cfg.polarity .* voltage);
  % prevent writing the same voltage over and over again
  previous_voltage = voltage;
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function control_auto_null(fig)
ft_info('auto null')
cfg     = getappdata(fig, 'cfg');
voltage  = getappdata(fig, 'voltage');
field   = getappdata(fig, 'field');
calib   = getappdata(fig, 'calib');

if ~isempty(calib)
  % this is incompatible with pairwise linked coils
  set(findobj(fig, 'tag', 'c1x'), 'value', 0);
  set(findobj(fig, 'tag', 'c3x'), 'value', 0);
  set(findobj(fig, 'tag', 'c5x'), 'value', 0);
  cb_enable_gui(fig); % update the gui

  residual = field(1:3);  % only keep the x, y, and z component
  residual = residual(:); % it should be a column vector
  voltage  = voltage(:);  % it should be a column vector

  % The idea is to compute a correction to the currently applied voltage that will bring the measured field to zero
  %   field = residual - calib * voltage  % this is the true environmental field
  %   residual = field + calib * voltage  % this is the residual field that we measure
  %
  % The sum of the environmental field plus that due to the voltage and the correction should be zero
  %   zero = field + calib * voltage + calib * correction
  %        = residual                + calib * correction
  % Hence
  %   residual = - calib * correction

  correction = - pinv(calib) * residual;
  voltage = voltage + correction;

  % clip the voltage to the specified range
  if any(voltage<cfg.range(1) | voltage>cfg.range(2))
    ft_info('clipping voltage')
    voltage(voltage>cfg.range(2)) = cfg.range(2);
    voltage(voltage<cfg.range(1)) = cfg.range(1);
  end

  % output the updated voltage
  setappdata(fig, 'voltage', voltage);
  control_voltage(fig);

else
  warning('cannot auto null, calibration has not yet been performed')
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibration_start_signal(fig)
% output a sine-wave with a different frequency on each of the coils
cfg     = getappdata(fig, 'cfg');
coils   = getappdata(fig, 'coils');
voltage = getappdata(fig, 'voltage');
fsample = cfg.calibration.fsample;
nsample = fsample*1; % one second of data is enough, it will repeat automatically
ncoil   = 6;

assert(numel(voltage)==ncoil);
assert(numel(cfg.calibration.frequency)==ncoil);
assert(numel(cfg.calibration.amplitude)==1);

coils.Rate = fsample;

% create one second of a sine-wave signal, it will loop until finished
voltage = zeros(ncoil, nsample);
time = (0:(nsample-1))/fsample;
for i=1:ncoil
  voltage(i,:) = cfg.polarity(i) * (cfg.calibration.amplitude * sin(cfg.calibration.frequency(i)*2*pi*time) + voltage(i));
end

if ~isempty(coils) && ~coils.Running
  preload(coils, voltage');
  start(coils, 'repeatoutput');
  while ~coils.Running
    pause(0.05);
  end
  fprintf('started calibration signal for %.1f seconds\n', cfg.calibration.duration + 2*cfg.calibration.pad);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibration_stop_signal(fig)
coils = getappdata(fig, 'coils');
if ~isempty(coils) && coils.Running
  stop(coils);
  while coils.Running
    pause(0.05);
  end
  flush(coils);
  ft_info('stopped calibration signal');
end
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calibration_compute(fig, field)
cfg     = getappdata(fig, 'cfg');
fsample = cfg.fsample;     % the sampling rate of the measured field
nsample = size(field, 2);  % the number of samples of the measured field
ncoil   = 6;

time = (0:(nsample-1))/fsample;
voltage = zeros(ncoil*2, nsample);
for i=1:ncoil
  % the output can be out of phase with the input, hence a cosine component is added
  voltage(2*i-1,:) = cfg.polarity(i) * cfg.calibration.amplitude * sin(cfg.calibration.frequency(i)*2*pi*time);
  voltage(2*i  ,:) = cfg.polarity(i) * cfg.calibration.amplitude * cos(cfg.calibration.frequency(i)*2*pi*time);
end

% drop the padding at the start and end
padsmp = (cfg.calibration.pad*fsample);
sel     = (padsmp+1):(nsample-padsmp);
time    = time(sel);
voltage = voltage(:,sel);
field   = field(:,sel);

% remove the DC component from the measured field
field = ft_preproc_baselinecorrect(field);

% compute the linear mix from the input voltages on the coils towards the measured output fields on the fluxgate
% field = calib * voltage + residue
calib = field / voltage;

model = calib*voltage;
residue = field - model;
gof = 1 - norm(residue, 'fro')/norm(field, 'fro');

figure; 
subplot(3,1,1); plot(time, field);   title('field');
subplot(3,1,2); plot(time, model);   title('model');
subplot(3,1,3); plot(time, residue); title('residue');
drawnow

% the measurement can be out of phase with the calibration signal
% hence the model voltage contains both a sine and a cosine component
c_sin = calib(1:2:end,:);
c_cos = calib(2:2:end,:);

% compute the phase difference between input and output
phase = atan2(c_cos, c_sin) * 180/pi;

% from the phase we can compute the delay in seconds
delay = phase;
for i=1:ncoil
  delay = (1./cfg.calibration.frequency) .* phase(i,:) / 360;
end

% combine the in-phase sine and the out-of-phase cosine component, assume that the delay is small
calib = sqrt(c_sin.^2 + c_cos.^2) * sign(v_sin);

fprintf('calibration computed\n');

disp('------------------ calibration in nT/V ----------------------')
disp(calib*1e9);
disp('------------------ phase in deg ----------------------')
disp(phase);
disp('------------------ delay in ms ----------------------')
disp(delay * 1000);

fprintf('goodness of fit = %.02f %%\n', gof*100);

setappdata(fig, 'calib', calib);
cb_enable_gui(fig);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measure_sample(fluxgate, varargin)
persistent fig counter nsample nchan buffer

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
    calibration_start_signal(fig);
    cfg = getappdata(fig, 'cfg');
    % include padding at the start and end
    nsample = round((cfg.calibration.duration + 2*cfg.calibration.pad) * cfg.fsample);
    nchan = numel(dat);
    buffer = nan(nchan, nsample);
    counter = 0;
    setappdata(fig, 'calibration', 'on');
  case 'on'
    if counter<nsample
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
end % switch calibration

% closed-loop updating of the voltage
switch getappdata(fig, 'closedloop')
  case 'on'
    calib  = getappdata(fig, 'calib');
    voltage = getappdata(fig, 'voltage');
    if ~isempty(calib)
      residual = dat(1:3);    % only keep the x, y, and z component
      residual = residual(:); % it should be a column vector
      voltage  = voltage(:);  % it should be a column vector

      correction = - pinv(calib) * residual;
      voltage = voltage + correction;
      setappdata(fig, 'voltage', voltage);
      control_voltage(fig);
    end
  case 'off'
    % nothing to do
end % switch running

% add the fluxgate field strength to the figure
setappdata(fig, 'field', dat(2:5));

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function measure_cleanup(fluxgate)
ft_info('cleanup');
configureCallback(fluxgate, 'off');
ft_info('stopped');
clear fluxgate
clear readSerialData
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_create_gui(h, eventdata, handles)
fig = getparent(h);

p1 = uipanel(fig, 'units', 'normalized', 'position', [0.0 0.2 0.5 0.8]);
p2 = uipanel(fig, 'units', 'normalized', 'position', [0.5 0.2 0.5 0.8]);
p3 = uipanel(fig, 'units', 'normalized', 'position', [0.0 0.0 1.0 0.2]);

uicontrol('parent', p1, 'style', 'text', 'string', 'voltage', 'units', 'normalized', 'position', [0.01 0.89 0.5 0.1], 'horizontalalignment', 'left');
uicontrol('parent', p2, 'style', 'text', 'string', 'field',   'units', 'normalized', 'position', [0.01 0.89 0.5 0.1], 'horizontalalignment', 'left');

uicontrol('parent', p1, 'tag', 'c1l', 'style', 'text', 'string', 'x');
uicontrol('parent', p1, 'tag', 'c1e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c1m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c1p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c2e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c2m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c2p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c1x', 'style', 'checkbox'); % link pairwise coils

uicontrol('parent', p1, 'tag', 'c3l', 'style', 'text', 'string', 'y');
uicontrol('parent', p1, 'tag', 'c3e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c3m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c3p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c4e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c4m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c4p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c3x', 'style', 'checkbox'); % link pairwise coils

uicontrol('parent', p1, 'tag', 'c5l', 'style', 'text', 'string', 'z');
uicontrol('parent', p1, 'tag', 'c5e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c5m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c5p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c6e', 'style', 'edit');
uicontrol('parent', p1, 'tag', 'c6m', 'style', 'pushbutton', 'string', '-');
uicontrol('parent', p1, 'tag', 'c6p', 'style', 'pushbutton', 'string', '+');
uicontrol('parent', p1, 'tag', 'c5x', 'style', 'checkbox'); % link pairwise coils
uicontrol('parent', p1, 'tag', 'cad', 'style', 'text', 'string', ''); % this is a dummy placeholder
uicontrol('parent', p1, 'tag', 'cax', 'style', 'checkbox');
uicontrol('parent', p1, 'tag', 'cal', 'style', 'text', 'string', 'cont. auto null', 'HorizontalAlignment', 'left');

ft_uilayout(p1, 'style', 'edit',       'callback', @cb_click);
ft_uilayout(p1, 'style', 'checkbox',   'callback', @cb_click);
ft_uilayout(p1, 'style', 'pushbutton', 'callback', @cb_click);

ft_uilayout(p1, 'tag', 'c..', 'units', 'normalized', 'height', 0.15)
ft_uilayout(p1, 'style', 'edit', 'width', 0.25, 'backgroundcolor', [1 1 1]);
ft_uilayout(p1, 'tag', 'c.p', 'width', 0.05)
ft_uilayout(p1, 'tag', 'c.m', 'width', 0.05)
ft_uilayout(p1, 'tag', 'c.x', 'width', 0.05)

ft_uilayout(p1, 'tag', 'c[12].', 'hpos', 'auto', 'vpos', 0.7);
ft_uilayout(p1, 'tag', 'c[34].', 'hpos', 'auto', 'vpos', 0.5);
ft_uilayout(p1, 'tag', 'c[56].', 'hpos', 'auto', 'vpos', 0.3);
ft_uilayout(p1, 'tag', 'ca.',    'hpos', 'auto', 'vpos', 0.1);
ft_uilayout(p1, 'tag', 'cal',    'width', 0.5);

uicontrol('parent', p2, 'tag', 'f1l', 'style', 'text', 'string', 'x');
uicontrol('parent', p2, 'tag', 'f1e', 'style', 'edit');
uicontrol('parent', p2, 'tag', 'f2l', 'style', 'text', 'string', 'y');
uicontrol('parent', p2, 'tag', 'f2e', 'style', 'edit');
uicontrol('parent', p2, 'tag', 'f3l', 'style', 'text', 'string', 'z');
uicontrol('parent', p2, 'tag', 'f3e', 'style', 'edit');
uicontrol('parent', p2, 'tag', 'f4l', 'style', 'text', 'string', 'abs');
uicontrol('parent', p2, 'tag', 'f4e', 'style', 'edit');

ft_uilayout(p2, 'tag', 'f..', 'units', 'normalized', 'height', 0.15)
ft_uilayout(p2, 'style', 'edit', 'width', 0.5, 'backgroundcolor', [0.9 0.9 0.9]);
ft_uilayout(p2, 'tag', 'f1.', 'hpos', 'auto', 'vpos', 0.7);
ft_uilayout(p2, 'tag', 'f2.', 'hpos', 'auto', 'vpos', 0.5);
ft_uilayout(p2, 'tag', 'f3.', 'hpos', 'auto', 'vpos', 0.3);
ft_uilayout(p2, 'tag', 'f4.', 'hpos', 'auto', 'vpos', 0.1);

uicontrol('parent', p3, 'tag', 'bc', 'style', 'pushbutton', 'string', 'calibrate', 'callback', @cb_click);
uicontrol('parent', p3, 'tag', 'ba', 'style', 'pushbutton', 'string', 'auto null', 'callback', @cb_click);
uicontrol('parent', p3, 'tag', 'bo', 'style', 'pushbutton', 'string', 'coils off', 'callback', @cb_click);
uicontrol('parent', p3, 'tag', 'bq', 'style', 'pushbutton', 'string', 'quit',      'callback', @cb_quit);

ft_uilayout(p3, 'style', 'pushbutton', 'units', 'normalized', 'width', 0.35, 'height', 0.6, 'hpos', 'auto', 'vpos', 0.2);

cb_enable_gui(fig);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_enable_gui(h, eventdata)
% ft_info('enable gui');
fig = getparent(h);
set(findall(fig, 'type', 'uicontrol', 'style', 'text'),       'Enable', 'on');
set(findall(fig, 'type', 'uicontrol', 'style', 'edit'),       'Enable', 'on');
set(findall(fig, 'type', 'uicontrol', 'style', 'pushbutton'), 'Enable', 'on');
set(findall(fig, 'type', 'uicontrol', 'style', 'checkbox'),   'Enable', 'on');

if isempty(getappdata(fig, 'calib'))
  % auto nulling is not possible
  set(findall(fig, 'tag', 'ba'), 'Enable', 'off');
  set(findall(fig, 'tag', 'cax'), 'Enable', 'off');
  set(findall(fig, 'tag', 'cal'), 'Enable', 'off');
end

% during continuous auto nulling
h = findall(fig, 'tag', 'cax');
if h.Value
  % the coils should not be pairwise linked
  set(findobj(fig, 'tag', 'c1x'), 'value', 0);
  set(findobj(fig, 'tag', 'c3x'), 'value', 0);
  set(findobj(fig, 'tag', 'c5x'), 'value', 0);
  % these should be disabled
  set(findall(fig, 'tag', 'ba'), 'Enable', 'off');
  set(findall(fig, 'tag', 'bc'), 'Enable', 'off');
  set(findall(fig, 'tag', 'bo'), 'Enable', 'off');
end

% disable the pairwise second coil when linked
h = findall(fig, 'tag', 'c1x');
if h.Value
  ft_uilayout(fig, 'tag', 'c2.', 'enable', 'off')
else
  ft_uilayout(fig, 'tag', 'c2.', 'enable', 'on')
end

% disable the pairwise second coil when linked
h = findall(fig, 'tag', 'c3x');
if h.Value
  ft_uilayout(fig, 'tag', 'c4.', 'enable', 'off')
else
  ft_uilayout(fig, 'tag', 'c4.', 'enable', 'on')
end

% disable the pairwise second coil when linked
h = findall(fig, 'tag', 'c5x');
if h.Value
  ft_uilayout(fig, 'tag', 'c6.', 'enable', 'off')
else
  ft_uilayout(fig, 'tag', 'c6.', 'enable', 'on')
end

% disable manual control during continuous auto nulling
h = findall(fig, 'tag', 'cax');
if h.Value
  ft_uilayout(fig, 'tag', 'c1.', 'enable', 'off')
  ft_uilayout(fig, 'tag', 'c2.', 'enable', 'off')
  ft_uilayout(fig, 'tag', 'c3.', 'enable', 'off')
  ft_uilayout(fig, 'tag', 'c4.', 'enable', 'off')
  ft_uilayout(fig, 'tag', 'c5.', 'enable', 'off')
  ft_uilayout(fig, 'tag', 'c6.', 'enable', 'off')
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_disable_gui(h, eventdata)
% ft_info('disable gui');
fig = getparent(h);
set(findall(fig, 'type', 'uicontrol', 'style', 'edit'),       'Enable', 'off');
set(findall(fig, 'type', 'uicontrol', 'style', 'pushbutton'), 'Enable', 'off');
set(findall(fig, 'type', 'uicontrol', 'style', 'checkbox'),   'Enable', 'off');
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata, timerarg)
% disp('redraw');
if nargin==3
  if isvalid(timerarg)
    fig = timerarg;
  else
    % this can happen when closing the application
    return
  end
else
  fig = getparent(h);
end
voltage = getappdata(fig, 'voltage');
field   = getappdata(fig, 'field');

% deal with the pairwise linkage between coils
if get(findall(fig, 'tag', 'c1x'), 'value')
  voltage(2) = voltage(1);
end
if get(findall(fig, 'tag', 'c3x'), 'value')
  voltage(4) = voltage(3);
end
if get(findall(fig, 'tag', 'c5x'), 'value')
  voltage(6) = voltage(5);
end
setappdata(fig, 'voltage', voltage);

% update the output voltage
set(findall(fig, 'tag', 'c1e'), 'string', voltage(1));
set(findall(fig, 'tag', 'c2e'), 'string', voltage(2));
set(findall(fig, 'tag', 'c3e'), 'string', voltage(3));
set(findall(fig, 'tag', 'c4e'), 'string', voltage(4));
set(findall(fig, 'tag', 'c5e'), 'string', voltage(5));
set(findall(fig, 'tag', 'c6e'), 'string', voltage(6));

% update the fluxgate field
set(findall(fig, 'tag', 'f1e'), 'string', sprintf('%.03f nT', 1e9*field(1)));
set(findall(fig, 'tag', 'f2e'), 'string', sprintf('%.03f nT', 1e9*field(2)));
set(findall(fig, 'tag', 'f3e'), 'string', sprintf('%.03f nT', 1e9*field(3)));
set(findall(fig, 'tag', 'f4e'), 'string', sprintf('%.03f nT', 1e9*field(4)));

if strcmp(getappdata(fig, 'closedloop'), 'off')
  % update the output voltage
  control_voltage(fig);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_click(h, eventdata)
fig     = getparent(h);
cfg     = getappdata(fig, 'cfg');
voltage = getappdata(fig, 'voltage');

% make steps of 1.0 for the full output range of -10V to +10V
% make steps of 0.1 for an output range of -1V to +1V
step = (cfg.range(2)-cfg.range(1))/20;

switch fig.SelectionType
  case 'alt'       % ctrl-click, this does not work on macOS
    step = step * 0.01;
  case 'extend'    % shift-click
    step = step * 0.1;
  otherwise
    step = step * 1;
end

switch h.Tag
  case 'c1e'
    voltage(1) = str2double(h.String);
  case 'c2e'
    voltage(2) = str2double(h.String);
  case 'c3e'
    voltage(3) = str2double(h.String);
  case 'c4e'
    voltage(4) = str2double(h.String);
  case 'c5e'
    voltage(5) = str2double(h.String);
  case 'c6e'
    voltage(6) = str2double(h.String);
  case 'c1p'
    voltage(1) = voltage(1) + step;
  case 'c1m'
    voltage(1) = voltage(1) - step;
  case 'c2p'
    voltage(2) = voltage(2) + step;
  case 'c2m'
    voltage(2) = voltage(2) - step;
  case 'c3p'
    voltage(3) = voltage(3) + step;
  case 'c3m'
    voltage(3) = voltage(3) - step;
  case 'c4p'
    voltage(4) = voltage(4) + step;
  case 'c4m'
    voltage(4) = voltage(4) - step;
  case 'c5p'
    voltage(5) = voltage(5) + step;
  case 'c5m'
    voltage(5) = voltage(5) - step;
  case 'c6p'
    voltage(6) = voltage(6) + step;
  case 'c6m'
    voltage(6) = voltage(6) - step;
  case 'c1x'
    cb_enable_gui(fig); % disable/enable the pairwise linked coil
  case 'c3x'
    cb_enable_gui(fig); % disable/enable the pairwise linked coil
  case 'c5x'
    cb_enable_gui(fig); % disable/enable the pairwise linked coil
  case 'cax'
    if h.Value
      ft_info('closedloop on');
      setappdata(fig, 'closedloop', 'on')
    else
      ft_info('closedloop off');
      setappdata(fig, 'closedloop', 'off')
    end
    cb_enable_gui(fig); % disable/enable the manual control of coils

  case 'bc'
    ft_info('initiating calibration')
    setappdata(fig', 'calibration', 'init');
    cb_disable_gui(fig);
    % the calibration is further handled in the measure_sample callback
  case 'ba'
    control_auto_null(fig);
    voltage = getappdata(fig, 'voltage');
  case 'bo'
    ft_info('switching coils off');
    voltage(:) = 0;
    % disable continuous auto nulling
    set(findobj(fig, 'tag', 'cax'), 'value', 0);
    cb_enable_gui(fig);

end % switch

% set near-zero values to zero
voltage(abs(voltage)<10*eps) = 0;

% clip the voltage to the specified range
if any(voltage<cfg.range(1) | voltage>cfg.range(2))
  ft_info('clipping voltage')
  voltage(voltage>cfg.range(2)) = cfg.range(2);
  voltage(voltage<cfg.range(1)) = cfg.range(1);
end

% update the values in the GUI
setappdata(fig, 'voltage', voltage);
cb_redraw(fig);
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)
fig = getparent(h);

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = fig.userdata;
else
  % determine the key that was pressed on the keyboard
  key = parsekeyboardevent(eventdata);
end

% get focus back to figure
if ~strcmp(h.Type, 'figure')
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
ft_info('switching coils off');
voltage = getappdata(fig, 'voltage');
voltage(:) = 0;
setappdata(fig, 'voltage', voltage);
control_voltage(fig);
setappdata(fig, 'closedloop', 'off');
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
