function data = edf2fieldtrip(filename)

% EDF2FIELDTRIP reads data from a EDF file with channels that have a different
% sampling rates. It upsamples all data to the highest sampling rate and
% concatenates all channels into a raw data structure that is compatible with the
% output of FT_PREPROCESSING.
%
% Use as
%   data = edf2fieldtrip(filename)
%
% For reading EDF files in which all channels have the same sampling rate, you can
% use the standard procedure with FT_DEFINETRIAL and FT_PREPROCESSING.
%
% See also FT_PREPROCESSING, FT_DEFINETRIAL, FT_REDEFINETRIAL

% Copyright (C) 2015, Robert Oostenveld
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

hdr = ft_read_header(filename);
samplerate = unique(hdr.orig.SampleRate);

data = cell(size(samplerate));

for i=1:numel(samplerate)
  chanindx = find(hdr.orig.SampleRate==samplerate(i));
  fprintf('reading %d channels with %g Hz sampling rate\n', numel(chanindx), samplerate(i));

  % read the header and data for the selected channels
  hdr = ft_read_header(filename, 'chanindx', chanindx);
  dat = ft_read_data(filename, 'header', hdr);

  % construct a time axis, starting at 0 seconds
  time = ((1:(hdr.nTrials*hdr.nSamples)) - 1)./hdr.Fs;

  % make a raw data structure
  data{i}.hdr   = hdr;
  data{i}.label = hdr.label;
  data{i}.time  = {time}; % only single data segment
  data{i}.trial = {dat};  % only single data segment
end

[maxrate, maxindex] = max(samplerate);

% upsample the data to the highest sampling rate
for i=1:numel(samplerate)
  if i==maxindex
    continue
  end
  fprintf('upsampling %d channels from %g to %g Hz\n', numel(data{i}.label), samplerate(i), maxrate);

  cfg = [];
  cfg.time = data{maxindex}.time;
  data{i} = ft_resampledata(cfg, data{i});
end

% concatenate them into a single data structure
data = ft_appenddata(cfg, data{:});

% reorder the channels to the original order in the EDF file
origlabel     = cellstr(hdr.orig.Label);
[currentorder, origorder] = match_str(origlabel, data.label); % sorted according to the 1st input argument
data.label    = data.label(origorder);
data.trial{1} = data.trial{1}(origorder,:);

% annotate the manual operation in the data structure provenance
cfg = [];
cfg.comment = 'reordered the channels to the original order in the EDF file';
data = ft_annotate(cfg, data);

if isfield(data, 'hdr')
  % remove this, as otherwise it might be very confusing with the subselections
  data = rmfield(data, 'hdr');
end
