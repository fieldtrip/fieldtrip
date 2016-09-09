function [event] = read_nex_event(filename)

% READ_NEX_EVENT for Plexon *.nex file, supports NEX variable types:
%   marker, interval, and event
%
% Use as
%   [event] = read_nex_event(filename)
%
% The event.type used to select events in ft_trialfun_general is the
% variable name from the NEX file (hdr.varheader.name - not to be confused
% with hdr.varheader.type).
%
% The sample numbers returned in event.sample correspond with the
% timestamps, correcting for the difference in sampling frequency in the
% continuous LFP channels and the system sampling frequency. Assuming 40kHz
% sampling frequency for the system and 1kHz for the LFP channels, it is
%   event.sample = timestamp / (40000/1000);
% If there are no continuous variables in the file, the system sampling
% frequency is used throughout, so
%   event.sample = timestamp;
%
% See also READ_NEX_HEADER, READ_NEX_DATA

% Copyright (C) 2005-2007, Robert Oostenveld
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

hdr = read_nex_header(filename);
adindx = find(cell2mat({hdr.varheader.typ})==5);
if isempty(adindx)  % this would otherwise produce an error
  warning('No continuous variables found - using hdr.filheader.frequency');
  smpfrq = hdr.filheader.frequency;
else
  smpfrq = hdr.varheader(adindx(1)).wfrequency;
end

fid=fopen(filename,'r','ieee-le');
event = struct('sample',{},'value',{},'timestamp',{},'type',{}, ...
  'duration',{},'offset',{});

% find the channel with the strobed trigger ("marker")
mrkvarnum = find([hdr.varheader.typ] == 6);

if ~isempty(mrkvarnum)  % both of these cases would have produced errors
  if numel(mrkvarnum) > 1
    warning('file contains >1 "marker" variable, reading only the first');
    mrkvarnum = mrkvarnum(1);
  end
  status = fseek(fid,hdr.varheader(mrkvarnum).offset,'bof');
  if status < 0;  error('error with fseek');  end
  
  % read the time of the triggers
  dum = fread(fid,hdr.varheader(mrkvarnum).cnt,'int32');
  timestamp = dum;
  dum = dum ./(hdr.filheader.frequency./smpfrq);
  mrk.tim = round(dum);
  
  % read the value of the triggers
  status = fseek(fid,64,'cof');
  if status < 0;  error('error with fseek');  end
  dum = fread(fid, [hdr.varheader(mrkvarnum).mrklen, ...
    hdr.varheader(mrkvarnum).cnt], 'uchar');
  mrk.val = str2num(char(dum(1:5,:)')); %#ok<ST2NM> non-scalar
  
  % translate into an FCDC event structure
  Nevent = length(mrk.tim);
  for i=1:Nevent
    event(i).sample         = mrk.tim(i);
    event(i).value          = mrk.val(i);
    event(i).timestamp      = timestamp(i);
    event(i).type           = hdr.varheader(mrkvarnum).nam;
    event(i).duration       = 1;
    event(i).offset         = 0;
  end
end

% find interval channels
intvarnum = find([hdr.varheader.typ] == 2);

for int = 1:numel(intvarnum)
  status = fseek(fid,hdr.varheader(intvarnum(int)).offset,'bof');
  if status < 0;  error('error with fseek');  end

  % read the time of the triggers
  dum1 = fread(fid,hdr.varheader(intvarnum(int)).cnt,'int32');
  dum2 = fread(fid,hdr.varheader(intvarnum(int)).cnt,'int32');
  timestamp = dum1;
  dum1 = dum1 ./(hdr.filheader.frequency./smpfrq);
  dum2 = dum2 ./(hdr.filheader.frequency./smpfrq);
  evt.tim = round(dum1);
  evt.dur = round(dum2-dum1);
  
  % translate into an FCDC event structure
  Nold = length(event);
  Nevent = length(evt.tim);
  for i=1:Nevent
    event(Nold+i).sample        = evt.tim(i);
    event(Nold+i).value         = [];
    event(Nold+i).timestamp     = timestamp(i);
    event(Nold+i).type          = hdr.varheader(intvarnum(int)).nam;
    event(Nold+i).duration      = evt.dur(i);
    event(Nold+i).offset        = 0;
  end
end

% find event channels
evtvarnum = find([hdr.varheader.typ] == 1);

for ev = 1:numel(evtvarnum)
  status = fseek(fid,hdr.varheader(evtvarnum(ev)).offset,'bof');
  if status < 0;  error('error with fseek');  end

  % read the time of the triggers
  dum = fread(fid,hdr.varheader(evtvarnum(ev)).cnt,'int32');
  timestamp = dum;
  dum = dum ./(hdr.filheader.frequency./smpfrq);
  evt.tim = round(dum);
  
  % translate into an FCDC event structure
  Nold = length(event);
  Nevent = length(evt.tim);
  for i=1:Nevent
    event(Nold+i).sample        = evt.tim(i);
    event(Nold+i).value         = [];
    event(Nold+i).timestamp     = timestamp(i);
    event(Nold+i).type          = hdr.varheader(evtvarnum(ev)).nam;
    event(Nold+i).duration      = 1;
    event(Nold+i).offset        = 0;
  end
end

status = fclose(fid);
if status < 0;  error('error with fclose');  end
