function [data, rem] = decode_modeeg(raw)

% DECODE_MODEEG takes a piece of the Modular EEG (OpenEEG) data
% stream from the serial port or bluetooth and decodes it.
%
% Use as
%  [dat, rem] = decode_modeeg(raw)
% where
%   raw = vector with bytes that was read from the serial port
%   dat = 6xN matrix with EEG values
%   rem = vector with remaining bytes that were not decoded
%
% The remaining bytes should be added to the raw data upon the next
% call.

% Copyright (C) 2012, Robert Oostenveld
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

packetsize = 17;


% ensure that it is an uint8 vector
if ~isa(raw, 'uint8')
  raw = uint8(raw);
end

begbyte   = find(raw(1:17)==hex2dec('a5') & raw(2:18)==hex2dec('5a'));
numpacket = floor((length(raw)-begbyte)/packetsize);
endbyte   = begbyte + numpacket*packetsize -1;

% trim the junk at the beginning
rem = raw(endbyte+1:end);
raw = raw(begbyte:endbyte);
raw = reshape(raw, 17, numpacket);

%   uint8_t   sync0;     // = 0xA5
%   uint8_t   sync1;     // = 0x5A
%   uint8_t   version;   // = 2
%   uint8_t   count;     // packet counter. Increases by 1 each packet
%   uint16_t  data[6];   // 10-bit sample (= 0 - 1023) in big endian (Motorola) format
%   uint8_t   switches;  // State of PD5 to PD2, in bits 3 to 0

sync0    = raw(1,:);
sync1    = raw(2,:);
version  = raw(3,:);
count    = raw(4,:);
data     = raw(5:16,:);
switches = raw(17,:);

% convert the data from uint8 into uint16 and swap the bytes
% this assumes that MATLAB is running on a little endian computer
data = swapbytes(reshape(typecast(data(:), 'uint16'), 6, numpacket));

