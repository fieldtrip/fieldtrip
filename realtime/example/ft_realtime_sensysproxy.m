function ft_realtime_sensysproxy(cfg)

% FT_REALTIME_SENSYSPROXY reads real-time continuous fluxgate magnetometer
% data from the serial interface provided by the Sensis FGM3D TD
% application and writes the data to the FieldTrip buffer.
%
% Use as
%   ft_realtime_sensysproxy(cfg)
%
% The configuration should contain
%   cfg.serialport  = string, name of the serial port (default = 'COM2')
%   cfg.blocksize   = number, in seconds (default = 0.1)
%   cfg.fsample     = sampling frequency (default = 1000)
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

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

cfg = ft_checkconfig(cfg);

% set the defaults
cfg.serialport        = ft_getopt(cfg, 'serialport', 'COM2');
cfg.baudrate          = ft_getopt(cfg, 'baudrate', 921600);
cfg.blocksize         = ft_getopt(cfg, 'blocksize', 0.1);
cfg.fsample           = ft_getopt(cfg, 'fsample', 400);
cfg.target            = ft_getopt(cfg, 'target');
cfg.target.datafile   = ft_getopt(cfg.target, 'datafile', 'buffer://localhost:1972');
cfg.target.dataformat = ft_getopt(cfg.target, 'dataformat', []);

% these are not used at the moment, the defaults seem to work fine
cfg.databits    = 8;
cfg.flowcontrol = 'none';
cfg.stopbits    = 1;
cfg.parity      = 'none';

% share the configuration details with the callback function
setappdata(0, 'cfg', cfg)

% the Sensys FGM3D TD Application allows serial data to be written out to a COM port
% open the serial device
sensys = serialport(cfg.serialport, cfg.baudrate);

%%
% this takes care of cleanup

cleanup = onCleanup(@()myCleanupFun(sensys));

% make sure the persistent variables are not reused from the last call
clear readSerialData

%%
% start the callback function, it processes each line

configureTerminator(sensys, "LF");
configureCallback(sensys, "terminator", @readSerialData);

% keep running forever, or until Ctrl-C
stopwatch = tic;
while (toc(stopwatch)<Inf) && strcmp(sensys.Status, 'open')
    pause(0.1);
end

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function readSerialData(sensys, event)

% this is called on every line/sample that is received
% but the data is not written that fast, so we need a buffer
persistent buffer nchan blocksize counter

cfg = getappdata(0, 'cfg');

if isempty(buffer)
    nchan = 5;
    % ensure that each block is an integer number of samples
    blocksize = round(cfg.blocksize*cfg.fsample);
    buffer = zeros(nchan, blocksize);
    counter = 0;
end

% The format is as follows:
% 1353411242214;0.0010416524;-0.0010926605;0.0014391395;0.0020856819<CR><LF>
%
% The separate values have the following meanings:
% Timestamp;Value Channel 1;Value Channel 2;Value Channel 3;Absolute Value<CR><LF>

line = char(readline(sensys));
line = strrep(line, ',', '.'); % depending on the language settings there might be a . or a , as decimal separator
dat = str2double(split(line, ';'));

buffer(:,rem(counter,blocksize)+1) = dat';
counter = counter + 1;

if counter==blocksize
    % after the first block
    hdr = [];
    hdr.Fs = cfg.fsample;
    hdr.nChans = nchan;
    hdr.nSamples = counter;
    hdr.nSamplesPre = 0;
    hdr.label = {'timestamp', 'x', 'y', 'z', 'abs'};
    % flush the file, write the header and subsequently write the data segment
    ft_write_data(cfg.target.datafile, buffer, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
    fprintf('wrote header and %d channels, %d samples\n', nchan, counter);
    
elseif mod(counter,blocksize)==0
    % after each subsequent block
    % write the data segment
    ft_write_data(cfg.target.datafile, buffer, 'append', true);
    fprintf('wrote %d channels, %d samples\n', nchan, counter);

end % process complete block

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myCleanupFun(sensys)
disp('cleanup');
configureCallback(sensys, "off");
disp('stopped')
clear sensys
clear readSerialData
end % function
