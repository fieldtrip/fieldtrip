function varargout = sccn_xdf(filename, hdr, begsample, endsample, chanindx)

% This is a wrapper to the reading function from the XDF MATLAB toolbox.
%
% Use as
%   hdr = sccn_xdf(filename);
%   dat = sccn_xdf(filename, hdr, begsample, endsample, chanindx);
%   evt = sccn_xdf(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, XDF2FIELDTRIP

% Copyright (C) 2019, Robert Oostenveld
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

% ensure this is on the path
ft_hastoolbox('xdf', 1);

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

streams = load_xdf(filename);

iscontinuous = false(size(streams));

% figure out which streams contain continuous/regular and discrete/irregular data
for i=1:numel(streams)
    
    % if the nominal srate is non-zero, the stream is considered continuous
    if ~strcmpi(streams{i}.info.nominal_srate, '0')
        
       iscontinuous(i) =  true;  
       
       if ~isfield(streams{i}.info, 'effective_srate') 
          
           % in case effective srate field value is missing, add one
           num_samples  = numel(streams{i}.time_stamps);
           t_begin      = streams{i}.time_stamps(1);
           t_end        = streams{i}.time_stamps(end);
           duration     = t_end - t_begin;
           streams{i}.info.effective_srate = (num_samples - 1) / duration;
           
           
       elseif isempty(streams{i}.info.effective_srate)
           
           % in case effective srate field value is missing, add one
           num_samples  = numel(streams{i}.time_stamps);
           t_begin      = streams{i}.time_stamps(1);
           t_end        = streams{i}.time_stamps(end);
           duration     = t_end - t_begin;
           streams{i}.info.effective_srate = (num_samples - 1) / duration;
           
       end 
       
    end

end

% determine the stream with the highest sampling rate
srate = nan(size(streams));
for i=1:numel(streams)
  if iscontinuous(i)
    srate(i) = streams{i}.info.effective_srate;
  end
end
[dum, indx] = max(srate);

% only keep the stream with the maximum sampling rate
% this is probably the EEG stream
stream = streams{indx};

if needhdr
  % this section of code is shared with xdf2fieldtrip
  hdr = [];
  if isfield(stream.info, 'effective_srate')
    % the stream contains continuously sampled data
    hdr.Fs                  = stream.info.effective_srate;
    hdr.nSamplesPre         = 0;
    hdr.nSamples            = length(stream.time_stamps);
    hdr.nTrials             = 1;
    hdr.FirstTimeStamp      = stream.time_stamps(1);
    hdr.TimeStampPerSample  = (stream.time_stamps(end)-stream.time_stamps(1)) / (length(stream.time_stamps) - 1);
  else
    % the stream does not contain continuously sampled data
    hdr.Fs                  = NaN;
    hdr.nSamplesPre         = NaN;
    hdr.nSamples            = NaN;
    hdr.nTrials             = NaN;
    hdr.FirstTimeStamp      = NaN;
    hdr.TimeStampPerSample  = NaN;
  end
  if isfield(stream.info.desc, 'channels')
    hdr.nChans    = numel(stream.info.desc.channels.channel);
  else
    hdr.nChans    = str2double(stream.info.channel_count);
  end
  hdr.label       = cell(hdr.nChans, 1);
  hdr.chantype    = cell(hdr.nChans, 1);
  hdr.chanunit    = cell(hdr.nChans, 1);
  
  prefix = stream.info.name;
  for j=1:hdr.nChans
    if isfield(stream.info.desc, 'channels')
      hdr.label{j} = [prefix '_' stream.info.desc.channels.channel{j}.label];
      hdr.chantype{j} = stream.info.desc.channels.channel{j}.type;
      hdr.chanunit{j} = stream.info.desc.channels.channel{j}.unit;
    else
      % the stream does not contain continuously sampled data
      hdr.label{j} = num2str(j);
      hdr.chantype{j} = 'unknown';
      hdr.chanunit{j} = 'unknown';
    end
  end
  
  % keep the original header details
  hdr.orig = stream.info;
  
  % return the header
  varargout = {hdr};
  
elseif needevt
  streams = streams(~iscontinuous);
  event = [];
  for i=1:numel(streams)
    for j=1:numel(streams{i}.time_series)
      
      % convert the timestamps to the corresponding sample in the selected data stream
      % the first sample in the data stream corresponds to 1
      timestamp = streams{i}.time_stamps(j);
      sample = round((timestamp - hdr.FirstTimeStamp)/hdr.TimeStampPerSample) + 1;
      
      event(end+1).type      = streams{i}.info.type;
      event(end  ).value     = streams{i}.time_series{j};
      event(end  ).timestamp = timestamp;
      event(end  ).sample    = sample;
    end
  end
  
  % return the events
  varargout = {event};
  
elseif needdat
  dat = stream.time_series(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
end
