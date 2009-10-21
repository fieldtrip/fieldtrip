function [data] = appendspike(cfg, varargin);

% APPENDSPIKE combines continuous data (i.e. LFP) with point-process data
% (i.e. spikes) into a single large dataset. For each spike channel an
% additional continuuos channel is inserted in the data that contains
% zeros most of the time, and an occasional one at the samples at which a
% spike occurred. The continuous and spike data are linked together using
% the timestamps.
%
% Use as
%   [spike] = appendspike(cfg, spike1, spike2, spike3, ...)
% where the input structures come from READ_FCDC_SPIKE, or as
%   [data]  = appendspike(cfg, data, spike1, spike2, ...)
% where the first data structure is the result of PREPROCESSING
% and the subsequent ones come from READ_FCDC_SPIKE.
%
% See also APPENDDATA, PREPROCESSING

% Copyright (C) 2007, Robert Osotenveld
%
% $Log: appendspike.m,v $
% Revision 1.4  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.3  2008/04/25 12:30:22  roboos
% enabled the checkdata call at the beginning
%
% Revision 1.2  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.1  2007/03/14 13:57:04  roboos
% initial implementation, tested on neuralynx data
%

fieldtripdefs

isspike = zeros(size(varargin));
for i=1:length(varargin)
  % this is a quick test, more rigourous checking is done later
  isspike(i) = isfield(varargin{i}, 'timestamp') & isfield(varargin{i}, 'label');
end

if all(isspike)
  spike = {};
  for i=1:length(varargin)
    spike{i} = checkdata(varargin{i}, 'datatype', 'spike');
  end

  % check the validity of the channel labels
  label = {};
  for i=1:length(spike)
    label = cat(1, label, spike{i}.label(:));
  end
  if length(unique(label))~=length(label)
    error('not all channel labels are unique');
  end

  % concatenate the spikes
  data = spike{1};
  for i=2:length(spike)
    data.label     = cat(2, data.label, spike{i}.label);
    data.waveform  = cat(2, data.waveform, spike{i}.waveform);
    data.timestamp = cat(2, data.timestamp, spike{i}.timestamp);
    data.unit      = cat(2, data.unit, spike{i}.unit);
  end

else
  % this checks the validity of the input data and simultaneously renames it for convenience
  data  = varargin{1}; % checkdata(varargin{1}, 'datatype', 'raw');
  spike = appendspike([], varargin{2:end}); 

  % check the validity of the channel labels
  label = cat(1, data.label(:), spike.label(:));
  if length(unique(label))~=length(label)
    error('not all channel labels are unique');
  end

  trl = findcfg(data.cfg, 'trl');
  if isempty(trl);
    error('could not find the trial information in the continuous data');
  end

  try
    FirstTimeStamp     = data.hdr.FirstTimeStamp;
    TimeStampPerSample = data.hdr.TimeStampPerSample;
  catch
    error('could not find the timestamp information in the continuous data');
  end

  for i=1:length(spike.label)
    % append the data with an empty channel
    data.label{end+1} = spike.label{i};
    for j=1:size(trl,1)
      data.trial{j}(end+1,:) = 0;
    end

    % determine the corresponding sample numbers for each timestamp
    ts = spike.timestamp{i};
    % timestamps can be uint64, hence explicitely convert to double at the right moment
    sample = round(double(ts-FirstTimeStamp)/TimeStampPerSample + 1);

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

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: appendspike.m,v 1.4 2008/09/22 20:17:43 roboos Exp $';
% remember the configuration details of the input data
cfg.previous = [];
for i=1:length(varargin)
  try, cfg.previous{i} = varargin{i}.cfg; end
end
% remember the exact configuration details in the output 
data.cfg = cfg;


