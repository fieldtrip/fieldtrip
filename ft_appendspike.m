function [data] = ft_appendspike(cfg, varargin)

% FT_APPENDSPIKE combines continuous data (i.e. LFP) with point-process data
% (i.e. spikes) into a single large dataset. For each spike channel an
% additional continuos channel is inserted in the data that contains
% zeros most of the time, and an occasional one at the samples at which a
% spike occurred. The continuous and spike data are linked together using
% the timestamps.
%
% Use as
%   [spike] = ft_appendspike(cfg, spike1, spike2, spike3, ...)
% where the input structures come from FT_READ_SPIKE, or as
%   [data]  = ft_appendspike(cfg, data, spike1, spike2, ...)
% where the first data structure is the result of FT_PREPROCESSING
% and the subsequent ones come from FT_READ_SPIKE.
%
% See also FT_APPENDDATA, FT_PREPROCESSING

% Copyright (C) 2007, Robert Osotenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance varargin

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

isspike = zeros(size(varargin));
for i=1:length(varargin)
  % this is a quick test, more rigourous checking is done later
  isspike(i) = isfield(varargin{i}, 'timestamp') & isfield(varargin{i}, 'label');
end

if all(isspike)
  spike = {};
  for i=1:length(varargin)
    % check if the input data is valid for this function
    spike{i} = ft_checkdata(varargin{i}, 'datatype', 'spike');
  end

  % check the validity of the channel labels
  label = {};
  for i=1:length(spike)
    label = cat(1, label, spike{i}.label(:));
  end
  if length(unique(label))~=length(label)
    ft_error('not all channel labels are unique');
  end

  % concatenate the spikes
  data = spike{1};
  for i=2:length(spike)
    data.label     = cat(2, data.label, spike{i}.label);

    % use a try construction in case a field is missing
    try, data.waveform  = cat(2, data.waveform, spike{i}.waveform); end
    try, data.timestamp = cat(2, data.timestamp, spike{i}.timestamp); end
    try, data.unit      = cat(2, data.unit, spike{i}.unit); end

    % these are optional fields, so use a try construction.
    try, data.time  = cat(2,data.time,spike{i}.time);    end
    try, data.trial = cat(2, data.trial,spike{i}.trial); end
    try, data.fourierspctrm = cat(2,data.fourierspctrm,spike{i}.fourierspctrm); end
  end

else
  % this checks the validity of the input data and simultaneously renames it for convenience
  data  = varargin{1}; % ft_checkdata(varargin{1}, 'datatype', 'raw');
  spike = ft_appendspike([], varargin{2:end});

  % check the validity of the channel labels
  label = cat(1, data.label(:), spike.label(:));
  if length(unique(label))~=length(label)
    ft_error('not all channel labels are unique');
  end

  if isfield(data, 'cfg')
    trl = ft_findcfg(data.cfg, 'trl');
  else
    trl = [];
  end

  if isempty(trl);
    ft_error('could not find the trial information in the continuous data');
  end

  try
    FirstTimeStamp     = data.hdr.FirstTimeStamp;
    TimeStampPerSample = data.hdr.TimeStampPerSample;
  catch
    ft_error('could not find the timestamp information in the continuous data');
  end

  for i=1:length(spike.label)
    % append the data with an empty channel
    data.label{end+1} = spike.label{i};
    for j=1:size(trl,1)
      data.trial{j}(end+1,:) = 0;
    end

    % determine the corresponding sample numbers for each timestamp
    ts = spike.timestamp{i};
    % timestamps can be uint64, hence explicitely convert to double at the
    % right moment
    if strcmp(class(ts),class(FirstTimeStamp))
      sample = round(double(ts-FirstTimeStamp)/TimeStampPerSample + 1);
    else
      sample = round(double(double(ts)-double(FirstTimeStamp))/TimeStampPerSample + 1);
    end
    fprintf('adding spike channel %s\n', spike.label{i});
    for j=1:size(trl,1)
      begsample = trl(j,1);
      endsample = trl(j,2);
      sel = find((sample>=begsample) & (sample<=endsample));
      fprintf('adding spike channel %s, trial %d contains %d spikes\n', spike.label{i}, j, length(sel));
      for k=1:length(sel)
        indx = sample(sel(k))-begsample+1;
        data.trial{j}(end,indx) = data.trial{j}(end,indx)+1;
      end % for k
    end % for j
  end % for i

end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous varargin
ft_postamble provenance data
ft_postamble history data
ft_postamble savevar data
