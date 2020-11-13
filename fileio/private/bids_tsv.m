function varargout = bids_tsv(filename, hdr, begsample, endsample, chanindx)

% BIDS_TSV reads time series data from a BIDS tsv and json file pair. This can for
% example be used to read the header and data from physio and stim files.
%
% Use as
%   hdr = bids_tsv(filename);
%   dat = bids_tsv(filename, hdr, begsample, endsample, chanindx);
%   evt = bids_tsv(filename, hdr);
% to read the header, the data or the event information.
%
% You should specify the name of the file containing the data as the filename, e.g.
% the _physio.tsv or the _stim.tsv file.
%
% See https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/06-physiological-and-other-continuous-recordings.html
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, MOTION_C3D, EVENTS_TSV

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

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

[p, f, x] = fileparts(filename);

datafile   = fullfile(p, [f '.tsv']);  % this can be a physio, stim, or another tsv file
headerfile = bids_sidecar(filename);   % this is the json file that describes the data
eventfile  = bids_sidecar(filename, 'events');

if needhdr || needdat
  fid = fopen_or_error(headerfile, 'r');
  json = fread(fid, [1, inf], 'char=>char');
  fclose(fid);
  
  % convert the json string into a structure
  json = jsondecode(json);
  
  % look at the first character of the file to determien whether it has a header line
  fid = fopen_or_error(datafile, 'r');
  str = fread(fid, 1, 'char=>char');
  fclose(fid);
  hasheader = ~isempty(regexp(str, '[azAZ]', 'once'));
  
  dat = readtable(datafile, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', hasheader);
  dat = table2array(dat)';
  assert(length(json.Columns)==size(dat,1), 'number of channels does not match');
end

if needhdr
  %% parse the header
  hdr.Fs = json.SamplingFrequency;
  hdr.nChans = length(json.Columns);
  hdr.nTrials = 1;
  hdr.nSamplesPre = 0;
  hdr.nSamples = size(dat,2);
  hdr.label = [json.Columns{:}];
  
  % remember the details from the json header file
  hdr.orig = json;
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% select the requested samples and channels
  dat = dat(chanindx, begsample:endsample);
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% read the events, this uses a function that can also be called directly from FT_READ_EVENT
  event = events_tsv(eventfile, hdr);
  
  % return the events
  varargout = {event};
  
end
