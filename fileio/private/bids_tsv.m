function varargout = bids_tsv(filename, hdr, begsample, endsample, chanindx)

% BIDS_TSV reads time series data from a BIDS tsv and json file pair
%
% Use as
%   hdr = bids_tsv(filename);
%   dat = bids_tsv(filename, hdr, begsample, endsample, chanindx);
%   evt = bids_tsv(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D

% Copyright (C) 2019 Robert Oostenveld
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
needdat = (nargin==5);

% use the full filename including path to distinguish between similarly named files in different directories
[p, f, x] = fileparts(filename);
if isempty(p)
  % no path was specified
  fullname = which(filename);
elseif startsWith(p, ['.' filesep])
  % a relative path was specified
  fullname = fullfile(pwd, p(3:end), [f, x]);
else
  fullname = filename;
end

[p, f, x] = fileparts(fullname);
jsonfile   = fullfile(p, [f '.json']);
tsvfile    = fullfile(p, [f '.tsv']);

if needhdr || needdat
  fid = fopen_or_error(jsonfile, 'r');
  str = fread(fid, [1, inf], 'char=>char');
  fclose(fid);
  % convert the json string into a structure
  orig = jsondecode(str);
  
  dat = readtable(tsvfile, 'FileType', 'text', 'Delimiter', '\t');
  dat = table2array(dat)';
  assert(length(orig.Columns)==size(dat,1), 'number of channels does not match');
end

if needhdr
  %% parse the header
  hdr.Fs = orig.SamplingFrequency;
  hdr.nChans = length(orig.Columns);
  hdr.nTrials = 1;
  hdr.nSamplesPre = 0;
  hdr.nSamples = size(dat,2);
  hdr.label = [orig.Columns{:}];
  
  % remember the details
  hdr.orig = orig;
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% parse the data
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% parse the events
  
  % construct the filename for the events
  [p, f, x] = fileparts(fullname);
  pieces = split(f, '_');
  base = sprintf('%s_', pieces{1:end-1}); % this ends with an inderscore
  eventsfile = fullfile(p, [base 'events.tsv']);
  
  if ~exist(eventsfile, 'file')
    ft_error('cannot read events');
  end
  
  % read the BIDS events as a table
  event = readtable(eventsfile, 'FileType', 'text', 'Delimiter', '\t');
  
  % Construct an events struct array that contains these fields
  %   event.type      = string
  %   event.sample    = expressed in samples, the first sample of a recording is 1
  %   event.value     = number or string
  %   event.offset    = expressed in samples
  %   event.duration  = expressed in samples
  %   event.timestamp = expressed in timestamp units, which vary over systems (optional)
  
  if any(strcmp(event.Properties.VariableNames, 'duration'))
    if iscell(event.duration)
      % it cannot have string values
      event.duration = [];
    elseif isnumeric(event.duration)
      % express the duration in samples
      %[event(:).duration] = deal(event.duration*hdr.Fs);
      event.duration = event.duration*hdr.Fs;
    end
  end
  
  event = table2struct(event);
  if ~isfield(event, 'sample')
    % the onset in events.tsv is always relative to the corresponding dataset
    onset = mat2cell([event(:).onset].*hdr.Fs, 1, ones(1,numel(event)));
    [event(:).sample] = deal(onset{:});
  end
  if ~isfield(event, 'type') && isfield(event, 'stim_type')
    type = {event(:).stim_type};
    [event(:).type] = deal(type{:});
  end
  timestamp = mat2cell([event(:).onset], 1, ones(1,numel(event)));
  [event(:).timestamp] = deal(timestamp{:});
  
  % remove all fields/columns that are not consistent with the output of FT_READ_EVENT
  event = keepfields(event, {'type', 'sample', 'value', 'offset', 'duration', 'timestamp'});
  
  % return the events
  varargout = {event};
  
end
