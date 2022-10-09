function ft_realtime_unicornproxy(cfg)

% FT_REALTIME_UNICORNPROXY reads continuous data from a Unicorn Hybrid Black
% wireless EEG system through Bluetooth and writes it to a FieldTrip buffer.
%
% The FieldTrip buffer is a network transparent server that allows the acquisition
% client to stream data to it. An analysis client can connect to read the data upon
% request. Multiple clients can connect simultaneously, each analyzing a specific
% aspect of the data concurrently.
%
% Use as
%   ft_realtime_unicornproxy(cfg)
%
% The configuration should contain
%   cfg.filename             = string, name of the serial port (default = '/dev/tty.FireFly-B106-SPP')
%   cfg.blocksize            = number, in seconds (default = 0.2)
%   cfg.writeeeg             = 'yes' or 'no' (default = 'yes')
%   cfg.writeaccel           = 'yes' or 'no' (default = 'no')
%   cfg.writegyro            = 'yes' or 'no' (default = 'no')
%   cfg.writebattery         = 'yes' or 'no' (default = 'no')
%   cfg.writecounter         = 'yes' or 'no' (default = 'no')
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C
%
% Prior to connecting the LED gives short flashes every second. After connecting it
% blinks on and off in a regular pace. When streaming the LED is constantly on.
%
% The Bluetooth protocol is documented in a PDF that is hosted on
% https://github.com/unicorn-bi/Unicorn-Suite-Hybrid-Black
%
% If you encounter Bluetooth connection problems on macOS, open a terminal and type
%   sudo pkill bluetoothd
%
% See also FT_REALTIME_SIGNALPROXY, FT_REALTIME_SIGNALVIEWER

% Copyright (C) 2022, Robert Oostenveld
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
cfg.filename          = ft_getopt(cfg, 'filename', '/dev/tty.UN-20211209');
cfg.baudrate          = ft_getopt(cfg, 'baudrate', 115200);
cfg.blocksize         = ft_getopt(cfg, 'blocksize', 0.2);
cfg.target            = ft_getopt(cfg, 'target');
cfg.target.datafile   = ft_getopt(cfg.target, 'datafile', 'buffer://localhost:1972');
cfg.target.dataformat = ft_getopt(cfg.target, 'dataformat', []);
cfg.writeeeg          = ft_getopt(cfg, 'writeeeg', 'yes');
cfg.writeaccel        = ft_getopt(cfg, 'writeaccel', 'no');
cfg.writegyro         = ft_getopt(cfg, 'writegyro', 'no');
cfg.writebattery      = ft_getopt(cfg, 'writebattery', 'no');
cfg.writecounter      = ft_getopt(cfg, 'writecounter', 'no');

% these are not used at the moment, the defaults seem to work fine
cfg.databits    = 8;
cfg.flowcontrol = 'none';
cfg.stopbits    = 1;
cfg.parity      = 'none';

% ensure that each block is an integer number of samples
cfg.blocksize = round(cfg.blocksize*250)/250;

% share the configuration details with the callback function
setappdata(0, 'cfg', cfg)

% the unicorn uses the SPP protocol, i.e. serial-over-bluetooth
% open the serial device
unicorn = serialport(cfg.filename, cfg.baudrate)


%%
% this takes care of cleanup

cleanup = onCleanup(@()myCleanupFun(unicorn));

%%

start_acq      = hex2dec({'0x61' '0x7C' '0x87'})';
start_response = hex2dec({'0x00' '0x00' '0x00'})';
stop_acq       = hex2dec({'0x63' '0x5C' '0xC5'})';
stop_response  = hex2dec({'0x00' '0x00' '0x00'})';
start_sequence = hex2dec({'0xC0' '0x00'})';
stop_sequence  = hex2dec({'0x0D' '0x0A'})';

%%
% start the data stream

write(unicorn, start_acq, 'uint8')

while unicorn.NumBytesAvailable<3
  % wait for the expected response
end
response = read(unicorn, 3, 'uint8');
assert(isequal(response, start_response));
disp('started')

%%
% start the callback function, it processes each block of 45 bytes
configureCallback(unicorn, "byte", 45, @readSerialData);

% keep running forever, or until Ctrl-C
stopwatch = tic;
while (toc(stopwatch)<Inf)
  pause(0.1);
end

end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function readSerialData(unicorn, event)

% this is called on every sample
% but the data is not written that fast, so we need a buffer
persistent buffer

cfg = getappdata(0, 'cfg');

% get a block of data
dat = uint8(read(unicorn, 45, 'uint8'));

% decipher the block of data
start = dat(1:2);
start_sequence = hex2dec({'0x0C0' '0x00'})';
assert(isequal(start, start_sequence))

% four bits
battery = bitand(dat(3), 0x0F);
battery = (100/15) * double(battery);

% big-endian, 24 bits
eeg = zeros(1,8);
for ch=1:8
  offset = (ch-1)*3;
  eeg1 = bitshift(uint32(dat(offset+4)), 16);
  eeg2 = bitshift(uint32(dat(offset+5)), 8);
  eeg3 = bitshift(uint32(dat(offset+6)), 0);
  %bitget(eeg1, 32:-1:1)
  %bitget(eeg2, 32:-1:1)
  %bitget(eeg3, 32:-1:1)
  eegv = bitor(bitor(eeg1, eeg2), eeg3);
  if bitget(eegv, 24)
    eegv = bitor(uint32(hex2dec('0xff000000')), eegv);
  end
  eegv = int32(eegv);
  %bitget(eeg, 32:-1:1)
  eeg(ch) = (4500000/50331642)*double(eegv);
end

% little-endian, 16 bits
accel = zeros(1,3);
for ch=1:3
  offset = (ch-1)*2;
  accel2 = bitshift(int16(dat(offset+28)), 0);
  accel1 = bitshift(int16(dat(offset+29)), 8);
  accelv = bitor(accel1, accel2);
  accel(ch) = (1/4096)*double(accelv);
end

% little-endian, 16 bits
gyro = zeros(1,3);
for ch=1:3
  offset = (ch-1)*2;
  gyro1 = bitshift(int16(dat(offset+34)), 0);
  gyro2 = bitshift(int16(dat(offset+35)), 8);
  gyrov = bitor(gyro1, gyro2);
  gyro(ch) = (1/32.8)*double(gyrov);
end

% little-endian, 32 bits
counter1 = bitshift(uint32(dat(40)), 0);
counter2 = bitshift(uint32(dat(41)), 8);
counter3 = bitshift(uint32(dat(42)), 16);
counter4 = bitshift(uint32(dat(43)), 24);
counter = double(bitor(counter1, bitor(counter2, bitor(counter3, counter4))));

% construct a representation that we can write to the FieldTrip buffer
dat = [];
hdr.Fs = 250;
hdr.label = {};

if istrue(cfg.writeeeg)
  dat = cat(2, dat, eeg);
  hdr.label = cat(2, hdr.label, {'EEG1' 'EEG2' 'EEG3' 'EEG4' 'EEG5' 'EEG6' 'EEG7' 'EEG8'});
end

if istrue(cfg.writeaccel)
  dat = cat(2, dat, accel);
  hdr.label = cat(2, hdr.label, {'AccelX' 'AccelY' 'AccelZ'});
end

if istrue(cfg.writegyro)
  dat = cat(2, dat, accel);
  hdr.label = cat(2, hdr.label, {'GyroX' 'GyroY' 'GyroZ'});
end

if istrue(cfg.writebattery)
  dat = cat(2, dat, battery);
  hdr.label = cat(2, hdr.label, {'Battery'});
end

if istrue(cfg.writecounter)
  dat = cat(2, dat, counter);
  hdr.label = cat(2, hdr.label, {'Counter'});
end

nchan = length(dat);
bufsize = 250*cfg.blocksize;

if counter==bufsize
  % after the first block
  % flush the file, write the header and subsequently write the data segment
  ft_write_data(cfg.target.datafile, buffer, 'header', hdr, 'dataformat', cfg.target.dataformat, 'append', false);
  fprintf('wrote %d channels, %d samples\n', nchan, counter);
elseif mod(counter,bufsize)==0
  % after each subsequent block
  % write the data segment
  ft_write_data(cfg.target.datafile, buffer, 'append', true);
  fprintf('wrote %d channels, %d samples\n', nchan, counter);
end % if count==1

% (re)initialize the buffer
if isempty(buffer) || mod(counter,bufsize)==0 || size(buffer,1)~=nchan
  buffer = zeros(nchan, bufsize);
end

% insert the data in the buffer
buffer(:,rem(counter,bufsize)+1) = dat';

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myCleanupFun(unicorn)
disp('cleanup');
configureCallback(unicorn, "off");
stop_acq = hex2dec({'0x63' '0x5C' '0xC5'})';
write(unicorn, stop_acq, 'uint8')
pause(0.1);
disp('stopped')
clear unicorn
clear readSerialData
end % function
