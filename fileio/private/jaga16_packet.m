function packet = jaga16_packet(buf, hastimestamp)

% JAGA16_PACKET converts the JAGA16 byte stream into packets

% Copyright (C) 2015-2016, Robert Oostenveld
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

buf = uint8(buf);
if buf(1)==0
  version = 0;
  buf = typecast(buf, 'uint16');
elseif buf(1)==3
  version = 3;
else
  error('unsupported version of packet');
end

% the packet size is 1396 bytes with timestamp, 1388 without.
if version==0 && hastimestamp
  % this is for the files created by the Jinga-Hi MATLAB and Python code
  buf = reshape(buf, 4+6+16*43, []);
  packet.ts1     = buf(1,:);
  packet.ts2     = buf(2,:);
  packet.ts3     = buf(3,:);
  packet.ts4     = buf(4,:);
  packet.ver     = buf(5,:);
  packet.nchan   = buf(6,:);
  packet.nbit    = buf(7,:);
  packet.fsample = buf(8,:);
  packet.sec     = buf(9,:);
  packet.smp     = buf(10,:);
  % the remainder is the data for this 16*43 block
  packet.dat = reshape(buf(11:end,:), 16, []);
elseif version==0
  % this is for the raw UDP data stream, where everything is shifted by four uint16 values (8 bytes)
  buf = reshape(buf, 0+6+16*43, []);
  packet.ver     = buf(1,:);
  packet.nchan   = buf(2,:);
  packet.nbit    = buf(3,:);
  packet.fsample = buf(4,:);
  packet.sec     = buf(5,:);
  packet.smp     = buf(6,:);
  % the remainder is the data for this 16*43 block
  packet.dat = reshape(buf(7:end,:), 16, []);
elseif version==3
  % this is for the raw UDP stream in version 3 format, each packet is 1388 bytes
  packet.ver    = buf(1);
  packet.nchan  = buf(2);
  packet.diagnostic_word  = typecast(buf([3 4]), 'uint16');
  packet.mode_word        = typecast(buf([5 6]), 'uint16');
  packet.fsample          = typecast(buf([7 8]), 'uint16');
  packet.smp              = typecast(buf([9:12]), 'uint32');
  % the remainder is the data for this 16*43 block
  packet.dat              = reshape(typecast(buf(13:end), 'uint16'), 16, []);
  % these are not present in v3
  packet.nbit   = 16;
  packet.sec    = NaN;
end

assert(all(packet.ver==0 | packet.ver==3));
assert(all(packet.nchan==16));
assert(all(packet.nbit==16));
