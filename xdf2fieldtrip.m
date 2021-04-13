function data = xdf2fieldtrip(filename, varargin)

% XDF2FIELDTRIP reads continuously sampled data from a XDF file with multiple
% streams. It upsamples the data of all streams to the highest sampling rate and
% concatenates all channels in all streams into a raw data structure that is
% compatible with the output of FT_PREPROCESSING.
%
% Use as
%   data = xdf2fieldtrip(filename, ...)
%
% Optional arguments should come in key-value pairs and can include
%   streamindx      = list, indices of the streams to read (default is all)
%   sraterange      = range of sampling rate in Hz in data streams to read [lowerbound, upperbound] 
%   streamkeywords  = cell array of keywords in stream names to index streams to read {keyword1, keyword2}
%
% You can also use the standard procedure with FT_DEFINETRIAL and FT_PREPROCESSING
% for XDF files. This will return (only) the continuously sampled stream with the
% highest sampling rate, which is typically the EEG.
%
% You can use FT_READ_EVENT to read the events from the non-continuous data streams.
% To get them aligned with the samples in one of the specific data streams, you
% should specify the corresponding header structure.
%
% See also FT_PREPROCESSING, FT_DEFINETRIAL, FT_REDEFINETRIAL

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

% process the options
streamindx      = ft_getopt(varargin, 'streamindx');
streamkeywords  = ft_getopt(varargin, 'streamkeywords');
sraterange      = ft_getopt(varargin, 'sraterange');

% for choosing streams, not more than one method can be selected
assert(sum(~cellfun(@isempty, {streamindx, streamkeywords, sraterange}))<=1, 'More than one method entered for stream selection: please use only one');

% ensure this is on the path
ft_hastoolbox('xdf', 1);

% read all streams
streams = load_xdf(filename);

% initialize an array of booleans indicating whether the streams are continuous
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

% give some feedback
for i=1:numel(streams)
  if iscontinuous(i)
    ft_info('stream %d contains continuous %s data\n', i, streams{i}.info.name);
  else
    ft_info('stream %d contains non-continuous %s data\n', i, streams{i}.info.name);
  end
end

% select the streams to continue working with
if isempty(streamindx) && isempty(streamkeywords)
  selected = true(size(streams));
else
  selected = false(size(streams));
  if ~isempty(streamindx)
      selected(streamindx) = true;
  end
  if ~isempty(streamkeywords)
      for si = 1:numel(streams)
          if contains(streams{si}.info.name, streamkeywords)
               selected(si) = true;
          end
      end
  end
end

% select the streams to continue working with by srate
if isempty(sraterange)
  inrange = true(size(streams));
else
  inrange = false(size(streams));  
  for i = 1:numel(streams)
      if isfield(streams{i}.info, 'effective_srate')
          if streams{i}.info.effective_srate >= sraterange(1) && streams{i}.info.effective_srate <= sraterange(2)
              inrange(i) = true;
          end
      end
  end
end


% discard the non-continuous streams
streams = streams(iscontinuous & selected & inrange);

if isempty(streams)
  ft_error('no continuous streams were selected');
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
  
  data{i}.hdr = hdr;
  data{i}.label = hdr.label;
  data{i}.time = {streams{i}.time_stamps};
  data{i}.trial = {streams{i}.time_series};
  
end % for all continuous streams

% determine the continuous stream with the highest sampling rate
srate = nan(size(streams));
for i=1:numel(streams)
  srate(i) = streams{i}.info.effective_srate;
end
[~, indx] = max(srate);


% copy the header from the stream with max srate
keephdr         = data{indx}.hdr;
keephdr.nChans  = 0; 
keephdr.label   = {};
keephdr.chantype = {};
keephdr.chanunit = {}; 

if numel(data)>1
    % resample all data structures, except the one with the max sampling rate
    % this will also align the time axes
    for i=1:numel(data)   
        
        % append channel information to the header
        keephdr.nChans      = keephdr.nChans + data{i}.hdr.nChans;
        keephdr.label       = [keephdr.label;       data{i}.hdr.label];
        keephdr.chantype    = [keephdr.chantype;    data{i}.hdr.chantype];
        keephdr.chanunit    = [keephdr.chanunit;    data{i}.hdr.chanunit];
                
        if i==indx
            continue
        end

        ft_notice('resampling %s', streams{i}.info.name);
        cfg = [];
        cfg.time = data{indx}.time;
        data{i} = ft_resampledata(cfg, data{i});
        
    end
    
    % append all data structures
    data = ft_appenddata([], data{:});
    
    % modify some fields in the header
    data.hdr = keephdr; 
else
  % simply return the first and only one
  data = data{1};
end
