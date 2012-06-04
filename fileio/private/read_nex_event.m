function [event] = read_nex_event(filename)

% READ_NEX_EVENT for Plexon *.nex file
%
% Use as
%   [event] = read_nex_event(filename)
%
% The sample numbers returned in event.sample correspond with the
% timestamps, correcting for the difference in sampling frequency in the
% continuous LFP channels and the system sampling frequency. Assuming 40kHz
% sampling frequency for the system and 1kHz for the LFP channels, it is
%   event.sample = timestamp / (40000/1000);
%
% See also READ_NEX_HEADER, READ_NEX_DATA

% Copyright (C) 2005-2007, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

hdr = read_nex_header(filename);
adindx = find(cell2mat({hdr.varheader.typ})==5);
smpfrq = hdr.varheader(adindx(1)).wfrequency;

% find the channel with the strobed trigger
mrkvarnum = find([hdr.varheader.typ] == 6);

fid=fopen(filename,'r','ieee-le');
status = fseek(fid,hdr.varheader(mrkvarnum).offset,'bof');

% read the time of the triggers
dum = fread(fid,hdr.varheader(mrkvarnum).cnt,'int32');
timestamp = dum;
dum = dum ./(hdr.filheader.frequency./smpfrq);
mrk.tim = round(dum);

% read the value of the triggers
status = fseek(fid,64,'cof');
dum = fread(fid,[hdr.varheader(mrkvarnum).mrklen,hdr.varheader(mrkvarnum).cnt],'uchar');
mrk.val = str2num(char(dum(1:5,:)'));

status = fclose(fid);

% translate into an FCDC event structure
Nevent = length(mrk.tim);
event = struct('sample', num2cell(mrk.tim), 'value', num2cell(mrk.val), 'timestamp', num2cell(timestamp));
for i=1:Nevent
  event(i).type     = hdr.varheader(mrkvarnum).nam;
  event(i).duration = 1;
  event(i).offset   = 0;
end

