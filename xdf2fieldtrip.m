function data = xdf2fieldtrip(filename)

% XDF2FIELDTRIP reads data from a XDF file with multiple streams. It upsamples the
% data of all streams to the highest sampling rate and concatenates all channels in
% all streams into a raw data structure that is compatible with the output of
% FT_PREPROCESSING.
%
% Use as
%   data = xdf2fieldtrip(filename)
%
% For reading XDF files with a single stream, you can use the standard procedure with
% FT_DEFINETRIAL and FT_PREPROCESSING.
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

% ensure this is on the path
ft_hastoolbox('xdf', 1);

streams = load_xdf(filename);

iscontinuous = false(size(streams));
% figure out which streams contain continuous/regular and discrete/irregular data
for i=1:numel(streams)
  iscontinuous(i) = isfield(streams{i}.info, 'effective_srate');
end

% discard the non-continuous streams
streams = streams(iscontinuous);

% convert each continuous stream into a FieldTrip raw data structure
data = cell(size(streams));
for i=1:numel(streams)
  if ischar(streams{i}.info.channel_count)
    streams{i}.info.channel_count = str2double(streams{i}.info.channel_count);
  end
  
  prefix = streams{i}.info.name;
  data{i}.label = cell(streams{i}.info.channel_count, 1);
  for j=1:streams{i}.info.channel_count
    data{i}.label{j} = [prefix '_' streams{i}.info.desc.channels.channel{j}.label];
  end
  data{i}.time = {streams{i}.time_stamps};
  data{i}.trial = {streams{i}.time_series};
end

% determine the continuous stream with the highest sampling rate
srate = nan(size(streams));
for i=1:numel(streams)
  srate(i) = streams{i}.info.effective_srate;
end
[~, indx] = max(srate);

% resample all data, except the one with the max sampling rate
for i=1:numel(data)
  if i==indx
    continue
  end
  cfg = [];
  cfg.time = data{indx}.time;
  data{i} = ft_resampledata(cfg, data{i});
end

% append all data
data = ft_appenddata([], data{:});
