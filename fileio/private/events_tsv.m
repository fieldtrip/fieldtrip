function varargout = events_tsv(filename, hdr)

% EVENTS_TSV is called from FT_READ_EVENT to read the events from a BIDS _events.tsv
% file. Although this function also reads the header for the sampling rate, it cannot
% be used to read data. Please see BIDS_TSV for that.
%
% Use as
%   hdr = events_tsv(filename)
%   evt = events_tsv(filename, hdr)
% to read the header or the event information.
%
% You should specify the _events.tsv file as the filename, the corresponding header
% file (with the sampling rate) will automatically be located in the same directory.
%
% See https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/05-task-events.html
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D, BIDS_TSV

% Copyright (C) 2018-2020, Robert Oostenveld
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

needhdr = (nargin==1);
needevt = (nargin==2);

if needhdr
  % find the data file and coresponding json file, these are needed for the sampling rate
  [datafile, jsonfile] = bids_datafile(filename);
  
  if isempty(jsonfile)
    ft_error('cannot find header information correspoding to ''%s''', filename);
  else
    ft_info('reading header information from ''%s''', jsonfile);
  end
  
  json = read_json(jsonfile);
  
  % The FieldTrip header has
  %   Fs, nSamples, nSamplesPre, nTrials, nChans, label, chantype, chanunit
  % and an optional json, elec, grad and opto field
  
  if isfield(json, 'SamplingFrequency')
    % this is for MEG, EEG and iEEG
    hdr.Fs = json.SamplingFrequency;
  elseif isfield(json, 'FrameRate')
    % this is for video, note that this can also include audio
    hdr.Fs = json.FrameRate;
  elseif isfield(json, 'AudioSampleRate')
    % this is for audio
    hdr.Fs = json.AudioSampleRate;
  else
    hdr.Fs = nan;
  end
  
  if isfield(json, 'RecordingDuration')
    % this is for MEG, EEG and iEEG
    hdr.nSamples = round(json.RecordingDuration*hdr.Fs);
  elseif isfield(json, 'Duration')
    % this is for video, note that this can also include audio
    hdr.nSamples = round(json.Duration*hdr.Fs);
  elseif isfield(json, 'AudioDuration')
    % this is for audio
    hdr.nSamples = round(json.AudioDuration*hdr.Fs);
  else
    hdr.nSamples = nan;
  end
  
  hdr.nSamplesPre  = 0;   % assume continuous
  hdr.nTrials      = 1;   % assume continuous
  hdr.nChans       = 0;   % number of channels
  hdr.label        = {};  % Nx1 cell-array with the label of each channel
  hdr.chantype     = {};  % Nx1 cell-array with the channel type, see FT_CHANTYPE
  hdr.chanunit     = {};  % Nx1 cell-array with the physical units, see FT_CHANUNIT
  
  fn = fieldnames(json);
  sel = find(endsWith(fn, 'ChannelCount'));
  for i=1:numel(sel)
    hdr.nChans = hdr.nChans + json.(fn{sel(i)});
  end
  for i=1:hdr.nChans
    hdr.label{i} = num2str(i);
    hdr.chantype{i} = 'unknown';
    hdr.chanunit{i} = 'unknown';
  end
  
  % keep the jsoninal header details
  hdr.json = json;
  
  % return the header
  varargout = {hdr};
end % if needhdr

if needevt
  opts = detectImportOptions(filename,'filetype','text', 'Delimiter', {'\t'} );
  % for some columns we have clear expectations on what they contain
  if any(strcmp(opts.VariableNames, 'onset'))
    opts = setvartype(opts, 'onset', 'double');
  end
  if any(strcmp('duration', opts.VariableNames))
    opts = setvartype(opts, 'duration', 'double');
  end
  if any(strcmp('sample', opts.VariableNames))
    opts = setvartype(opts, 'sample', 'double');
  end
  if any(strcmp('offset', opts.VariableNames))
    opts = setvartype(opts, 'offset', 'double');
  end
  if any(strcmp(opts.VariableNames, 'type'))
    opts = setvartype(opts, 'type', 'char');
  end
  if any(strcmp(opts.VariableNames, 'trial_type'))
    opts = setvartype(opts, 'trial_type', 'char');
  end
  if any(strcmp(opts.VariableNames, 'event_type'))
    opts = setvartype(opts, 'event_type', 'char');
  end
  if any(strcmp(opts.VariableNames, 'stim_type'))
    opts = setvartype(opts, 'event_type', 'char');
  end
  % initially keep the value as a string, when possible it will be converted to a number later
  if any(strcmp(opts.VariableNames, 'value'))
    opts = setvartype(opts, 'value', 'char');
  end
  if any(strcmp(opts.VariableNames, 'event_value'))
    opts = setvartype(opts, 'event_value', 'char');
  end
  if any(strcmp(opts.VariableNames, 'stim_value'))
    opts = setvartype(opts, 'event_value', 'char');
  end
  
  % this keeps the type and value column as string, and castst others into doubles right away
  tsv = readtable(filename, opts);
  
  % The FieldTrip event structure should have
  %    type, value, sample, offset, duration
  % and an optional timestamp field.
  
  % start with an empty structure
  event = struct();
  
  % the event type should be a string
  if iscolumn(tsv, 'type')
    for k = 1:size(tsv,1)
      event(k).type = tsv.type{k};
    end
  elseif iscolumn(tsv, 'trial_type')
    for k = 1:size(tsv,1)
      event(k).type = tsv.trial_type{k};
    end
  elseif iscolumn(tsv, 'event_type')
    for k = 1:size(tsv,1)
      event(k).type = tsv.event_type{k};
    end
  elseif iscolumn(tsv, 'stim_type')
    for k = 1:size(tsv,1)
      event(k).type = tsv.stim_type{k};
    end
  else
    % assign the type for all events as empty
    event(1).type = [];
  end
  
  % the event value can be a string or numeric
  % when possible, the conversion of string to numeric will be done in FT_READ_EVENT
  if iscolumn(tsv, 'value')
    for k = 1:size(tsv,1)
      event(k).value = tsv.value{k};
    end
  elseif iscolumn(tsv, 'event_value')
    for k = 1:size(tsv,1)
      event(k).value = tsv.event_value{k};
    end
  elseif iscolumn(tsv, 'stim_value')
    for k = 1:size(tsv,1)
      event(k).value = tsv.stim_value{k};
    end
  else
    % assign the value for all events as empty
    event(1).value = [];
  end
  
  if iscolumn(tsv, 'sample')
    % use the specified sample number, these are assumed to be one-offset
    for k = 1:size(tsv,1)
      event(k).sample = tsv.sample(k);
    end
  elseif iscolumn(tsv, 'onset') && nargin>1
    % construct the sample number from the onset in seconds
    for k = 1:size(tsv,1)
      event(k).sample = round(tsv.onset(k)*hdr.Fs) + 1;
    end
  else
    % we don't know the sampling rate, so cannot determine the sample number
    event(1).sample = [];
  end
  
  if iscolumn(tsv, 'offset')
    for k = 1:size(tsv,1)
      event(k).offset = tsv.offset(k);
    end
  else
    % assign the offset for all events as empty
    event(1).offset = [];
  end
  
  if iscolumn(tsv, 'duration')
    for k = 1:size(tsv,1)
      % the onset and duration in the BIDS events.tsv file are expressed in samples, not in seconds
      event(k).duration = tsv.duration(k)*hdr.Fs;
    end
  else
    % assign the duration for all events as empty
    event(1).duration = [];
  end
  
  if iscolumn(tsv, 'timestamp')
    % use the specified timestamp, which is un unknown units
    for k = 1:size(tsv,1)
      event(k).timestamp = tsv.timestamp(k);
    end
  else
    % use the onset in seconds as the timestamp
    for k = 1:size(tsv,1)
      event(k).timestamp = tsv.onset(k);
    end
  end
  
  % return the events
  varargout = {event};
end % if needevt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = iscolumn(tsv, name)
bool = ismember(name, tsv.Properties.VariableNames);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function json = read_json(filename)
ft_hastoolbox('jsonlab', 1);
json = loadjson(filename);
json = ft_struct2char(json); % convert strings into char-arrays


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = is_a_number(input)
pattern = '[0-9\.]';
bool    = length(regexp(input, pattern))==length(input);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = is_a_nan(input)
pattern = 'nan|n/a';
bool    = ~isempty(regexp(lower(input), pattern, 'once'));
