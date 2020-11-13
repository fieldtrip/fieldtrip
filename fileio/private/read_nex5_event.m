function [event] = read_nex5_event(filename)

% READ_NEX5_EVENT for Nex Technologies *.nex5 file, supports NEX5 variable types:
%   marker, interval, and event
%
% Use as
%   [event] = read_nex5_event(filename)
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
% See also READ_NEX5_HEADER, READ_NEX5

% Copyright (C) 2020, Robert Oostenveld, Alex Kirillov
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

hdr = read_nex5_header(filename);
adindx = find(cell2mat({hdr.VarHeader.Type})==5);
if isempty(adindx)  % this would otherwise produce an error
  ft_warning('No continuous variables found - using hdr.FileHeader.Frequency');
  samplingFreq = hdr.FileHeader.Frequency;
else
  samplingFreq = hdr.VarHeader(adindx(1)).WFrequency;
end

divisorToSampleNumber = hdr.FileHeader.Frequency./samplingFreq;

fid = fopen_or_error(filename,'r','ieee-le');
event = struct('sample',{},'value',{},'timestamp',{},'type',{}, ...
  'duration',{},'offset',{});

% find any channels with strobed triggers ("markers")
mrkvarnum = find([hdr.VarHeader.Type] == 6);

for mrkn = 1:numel(mrkvarnum)
  vh = hdr.VarHeader(mrkvarnum(mrkn));
  status = fseek(fid,vh.DataOffset,'bof');
  if status < 0;  ft_error('error with fseek');  end
  
  % read the times of the triggers
  ts = Nex5ReadTimestampsAsColumn(fid, vh);
  
  % skip importing this marker if empty
  name = deblank(vh.Name);
  if isempty(ts)
    ft_warning(['skipping marker ' name ' because no timestamps were found'])
    continue
  end
  
  % convert timestamp to sample number
  timestamps = ts;
  ts = ts./divisorToSampleNumber;
  roundTimes = round(ts);
  
  nevents = length(ts);
  % read the values of the markers
  for j=1:vh.NumberOfMarkerFields
    markerName = deblank(fread(fid, [1 64], 'uint8=>char'));
    if vh.MarkerDataType == 0
        % marker values are sored as strings
        markersStored = fread(fid, [vh.MarkerLength vh.Count], 'uchar');
        markerStrings = deblank(char(markersStored(1:vh.MarkerLength, :))');
        numValues = str2num(markerStrings);
        if length(numValues) > 0
          markerEventValues = num2cell(numValues);
        else
          % empty numValues means that we cannot convert markerStrings to numbers
          markerEventValues = cellstr(markerStrings);
        end
    else
      % marker values are sored as numbers
      markerEventValues = num2cell(fread(fid, vh.Count, 'uint32'));
    end
    tmp = struct('sample', num2cell(roundTimes), ...
      'value', markerEventValues, ...
      'timestamp', num2cell(timestamps), ...
      'type', repmat({[name '_' markerName]},[nevents,1]), ...
      'duration', num2cell(ones(nevents,1)), ...
      'offset', num2cell(zeros(nevents,1)));
    event = [event; tmp];
  end
end

% find interval channels
intvarnum = find([hdr.VarHeader.Type] == 2);

for intVarIndex = 1:numel(intvarnum)
  vh = hdr.VarHeader(intvarnum(intVarIndex));
  status = fseek(fid, vh.DataOffset,'bof');
  if status < 0;  ft_error('error with fseek');  end

  % read the time of the triggers
  intStarts = Nex5ReadTimestampsAsColumn(fid, vh);
  intEnds = Nex5ReadTimestampsAsColumn(fid, vh);
  
  % skip importing this interval if empty
  if isempty(intStarts) && isempty(intEnds)
    ft_warning(['skipping interval ' deblank(vh.Name) ' because no timestamps were found'])
    continue
  end
  
  timestamps = intStarts;
  intStarts = intStarts ./divisorToSampleNumber;
  intEnds = intEnds ./divisorToSampleNumber;
  intEvent.times = round(intStarts);
  intEvent.durations = round(intEnds - intStarts);
  
  % translate into an FCDC event structure
  nevents = length(intEvent.times);
  tmp = struct('sample',num2cell(intEvent.times), ...
    'value', cell(nevents,1), ...
    'timestamp', num2cell(timestamps), ...
    'type', repmat({vh.Name},[nevents,1]), ...
    'duration', num2cell(intEvent.durations), ...
    'offset', num2cell(zeros(nevents,1)));
  event = [event; tmp];
end

% find event channels
evtvarnum = find([hdr.VarHeader.Type] == 1);

for ev = 1:numel(evtvarnum)
  vh = hdr.VarHeader(evtvarnum(ev));
  status = fseek(fid, vh.DataOffset, 'bof');
  if status < 0;  ft_error('error with fseek');  end

  % read the times of the events
  ts = Nex5ReadTimestampsAsColumn(fid, vh);
  
  % skip importing this event if empty
  if isempty(ts)
    ft_warning(['skipping event ' deblank(vh.Name) ' because no timestamps were found'])
    continue
  end
  
  % convert timestamp to sample number
  timestamps = ts;
  ts = ts./divisorToSampleNumber;
  roundTimes = round(ts);
  
  % translate into an FCDC event structure
  nevents = length(ts);
  tmp = struct('sample', num2cell(roundTimes), ...
    'value', cell(nevents,1), ...
    'timestamp', num2cell(timestamps), ...
    'type', repmat({vh.Name},[nevents,1]), ...
    'duration', num2cell(ones(nevents,1)), ...
    'offset', num2cell(zeros(nevents,1)));
  event = [event; tmp];
end

status = fclose(fid);
if status < 0;  ft_error('error with fclose');  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction to read nex5 timestamps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ts = Nex5ReadTimestampsAsColumn(fid, varHeader)
  if varHeader.TimestampDataType == 0
    ts = fread(fid, varHeader.Count, 'int32');
  else
    ts = fread(fid, varHeader.Count, 'int64'); 
  end
end