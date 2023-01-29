function varargout = eegsynth_tsv(filename, hdr)

% EEGSYNTH_TSV is called from FT_READ_EVENT to read the events from a tsv file
% written by the recordtrigger module. The .tsv file should also contain a
% synchronization trigger from the recordsignal module.
%
% Use as
%   hdr = events_tsv(filename)
%   evt = events_tsv(filename, hdr)
% to read the header or the event information. Note that when reading the header, the
% number of channels in the actual data is unknown.
%
% See https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/05-task-events.html
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT, QUALISYS_TSV, EVENTS_TSV, BIDS_TSV

% FIXME it might be that the events are one sample off, but I cannot be bothered
% checking that precisely at the moment of the initial implementation.

% Copyright (C) 2022, Robert Oostenveld
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

% read the table from file
t = readtable(filename, 'filetype', 'text', 'delimiter', '\t');

type = unique(t.event);
sync = type(endsWith(type, 'synchronize'));
assert(length(sync)==1, 'there should be exactly one type of trigger that ends with "synchronize"')
sync = sync{1};

% find the mapping between timestamps and samples
sel = strcmp(t.event, sync);
firstTimeStamp = t.timestamp(1);
sec = seconds(t.timestamp(sel) - firstTimeStamp);
smp = t.value(sel); % the "synchronize" triggers contain the sample number of the corresponding EEG file

% figure
% plot(sec, smp, '.')
% xlabel('seconds')
% ylabel('samples')

if needhdr
  % fit a first-order polynomial that expresses samples per timestamp
  p = polyfit(sec, smp, 1);
  samplesPerTimeStamp = p(1);

  % fit a first-order polynomial that expresses timestamps per sample
  p = polyfit(smp, sec, 1);

  % construct a dummy header
  hdr = [];
  hdr.nChans = 1;
  hdr.label = {'unknown'};
  hdr.nSamples = inf;
  hdr.nSamplesPre = 0;
  hdr.nTrials = 1;
  hdr.Fs = samplesPerTimeStamp;
  hdr.FirstTimeStamp = firstTimeStamp + seconds(polyval(p, 1));
  hdr.TimeStampPerSample  = samplesPerTimeStamp;

  % return the header
  varargout = {hdr};
end % if needhdr

if needevt
  % fit a first-order polynomial that expresses samples per timestamp
  p = polyfit(sec, smp, 1);

  % estimate the sample number from the timestamp for all other events
  sec = seconds(t.timestamp - firstTimeStamp);
  smp = round(polyval(p, sec)) - 1;

  event = struct;
  for i=1:numel(smp)
    event(i).type = t.event{i};
    event(i).value = t.value(i);
    event(i).sample = smp(i);
    event(i).duration = [];
    event(i).offset = [];
  end

  % return the events
  varargout = {event};
end % if needevt
