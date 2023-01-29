function [data, event] = xdf2fieldtrip(filename, varargin)

% XDF2FIELDTRIP reads continuously sampled data from a XDF file with multiple
% streams. It upsamples the data of all streams to the highest sampling rate and
% concatenates all channels in all streams into a raw data structure that is
% compatible with the output of FT_PREPROCESSING.
%
% Use as
%   [data, events] = xdf2fieldtrip(filename, ...)
%
% Optional arguments should come in key-value pairs and can include
%   streamindx      = number or list, indices of the streams to read (default is all)
%   streamrate      = [lowerbound upperbound], read only data streams within this range of sampling rates (in Hz)
%   streamkeywords  = cell-array with strings, keywords contained in the stream to read
%
% You can also use the standard procedure with FT_DEFINETRIAL and FT_PREPROCESSING
% for XDF files. This will return (only) the continuously sampled stream with the
% highest sampling rate, which is typically the EEG.
%
% You can also use FT_READ_EVENT to read the events from the non-continuous data
% streams. To get them aligned with the samples in one of the specific data streams,
% you should specify the corresponding header structure.
%
% See also FT_PREPROCESSING, FT_DEFINETRIAL, FT_REDEFINETRIAL

% Copyright (C) 2019-2021, Robert Oostenveld
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

% process the options
streamindx      = ft_getopt(varargin, 'streamindx');
streamkeywords  = ft_getopt(varargin, 'streamkeywords');
streamrate      = ft_getopt(varargin, 'streamrate');

if isempty(streamrate)
  % this option used to be called sraterange
  % on 10 June 2021 it was changed into streamrate for consistency with the other options
  streamrate = ft_getopt(varargin, 'sraterange');
end

if length(streamrate)==1
  % it should be specified as [low high]
  streamrate = [streamrate streamrate];
end

% not more than one method can be selected for choosing streams
assert(sum(~cellfun(@isempty, {streamindx, streamkeywords, streamrate}))<=1, 'you can only use a single method for stream selection');

% ensure this is on the path
ft_hastoolbox('xdf', 1);

% read all streams
streams = load_xdf(filename);

% initialize an array of booleans indicating whether the streams are continuous
iscontinuous = false(size(streams));

% figure out which streams contain continuous/regular, and which ones contain discrete/irregular data
for i=1:numel(streams)
  % the stream is considered continuous if the nominal srate is non-zero
  if str2double(streams{i}.info.nominal_srate)~=0
    iscontinuous(i) = true;
    
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
    
  end % if nonzero nominal sampling rate
end % for all streams

% give some feedback
for i=1:numel(streams)
  if iscontinuous(i)
    ft_info('stream %d contains continuous %s data\n', i, streams{i}.info.name);
  else
    ft_info('stream %d contains non-continuous %s data\n', i, streams{i}.info.name);
  end
end

% select the streams to continue working based on keywords
if isempty(streamindx) && isempty(streamkeywords)
  haskeyword = true(size(streams));
else
  haskeyword = false(size(streams));
  if ~isempty(streamindx)
    haskeyword(streamindx) = true;
  end
  if ~isempty(streamkeywords)
    for i=1:numel(streams)
      haskeyword(i) = contains(streams{i}.info.name, streamkeywords);
    end
  end
end

% select the streams to continue working based on sampling rate
if isempty(streamrate)
  inrange = true(size(streams));
else
  inrange = false(size(streams));
  for i = 1:numel(streams)
    if isfield(streams{i}.info, 'effective_srate')
      if streams{i}.info.effective_srate >= streamrate(1) && streams{i}.info.effective_srate <= streamrate(2)
        inrange(i) = true;
      end
    end
  end
end

% convert the non-continuous streams to events
event = [];
for i=1:numel(streams)
  if iscontinuous(i)
    continue
  end
  for k=1:length(streams{i}.time_stamps)
    try
      event(end+1).type      = streams{1}.info.type;
      event(end  ).value     = streams{i}.time_series{k}; % this is a cell-array with strings
      event(end  ).sample    = nan;   % not defined, as it is not clear to which continuous stream the event relates
      event(end  ).duration  = [];    % not specified
      event(end  ).offset    = [];    % not specified
      event(end  ).timestamp = streams{i}.time_stamps(k); % this is a scalar array
    catch
      % skip this event, the formatting might be different than assumed in the code above
    end
  end
end

% continue with the continuous streams
streams = streams(iscontinuous & haskeyword & inrange);

if isempty(streams) && nargout==1
  % in case of two output arguments it will return the events
  ft_error('no continuous streams were present or selected');
end

% convert each continuous stream into a FieldTrip raw data structure
data = cell(size(streams));
for i=1:numel(streams)
  
  % make a copy for convenience
  stream = streams{i};
  
  % this section of code is shared with fileio/private/sccn_xdf
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
      if isfield(stream.info.desc.channels.channel{j}, 'label')
        hdr.label{j} = [prefix '_' stream.info.desc.channels.channel{j}.label];
      else
        hdr.label{j} = num2str(j);
      end
      if isfield(stream.info.desc.channels.channel{j}, 'type')
        hdr.chantype{j} = stream.info.desc.channels.channel{j}.type;
      else
        hdr.chantype{j} = 'unknown';
      end
      if isfield(stream.info.desc.channels.channel{j}, 'unit')
        hdr.chanunit{j} = stream.info.desc.channels.channel{j}.unit;
      else
        hdr.chanunit{j} = 'unknown';
      end
    else
      % the stream does not contain continuously sampled data
      hdr.label{j} = num2str(j);
      hdr.chantype{j} = 'unknown';
      hdr.chanunit{j} = 'unknown';
    end
  end
  
  % keep the original header details
  hdr.orig = stream.info;
  
  data{i}.hdr = hdr;
  data{i}.label = hdr.label;
  data{i}.time = {streams{i}.time_stamps};
  data{i}.trial = {streams{i}.time_series};
  
end % for all continuous streams


if numel(data)>1
  % determine the continuous stream with the highest sampling rate
  srate = nan(size(streams));
  for i=1:numel(streams)
    srate(i) = streams{i}.info.effective_srate;
  end
  [dum, highest] = max(srate);
  
  % copy the header from the stream with the highest sampling rate
  keephdr          = data{highest}.hdr;
  keephdr.nChans   = 0;
  keephdr.label    = {};
  keephdr.chantype = {};
  keephdr.chanunit = {};
  
  % resample all data structures, except the one with the max sampling rate
  % this will also align the time axes
  for i=1:numel(data)
    
    % append this stream channel information to the combined header
    keephdr.nChans   = keephdr.nChans + data{i}.hdr.nChans;
    keephdr.label    = [keephdr.label;       data{i}.hdr.label];
    keephdr.chantype = [keephdr.chantype;    data{i}.hdr.chantype];
    keephdr.chanunit = [keephdr.chanunit;    data{i}.hdr.chanunit];
    
    if i==highest
      continue
    end
    
    ft_notice('resampling %s', streams{i}.info.name);
    
    cfg = [];
    cfg.time = data{highest}.time;
    data{i} = ft_resampledata(cfg, data{i});
    
  end
  
  % append all data structures
  data = ft_appenddata([], data{:});
  
  % modify some fields in the header
  data.hdr = keephdr;
  
elseif numel(data)==1
  % simply return the first (and only) data structure
  data = data{1};
  
else
  % do not return any data, but possibly return events
  data = [];
end
