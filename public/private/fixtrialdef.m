function data = fixtrialdef(data)

% FIXTRIALDEF adds a trl matrix to the raw data configuration
% which is constructed on the fly, assuming that the trials are
% consecutive segments of a continuous recording

% Copyright (C) 2009-2010, Robert Oostenveld and Jan-Mathijs Schoffelen

if isfield(data, 'sampleinfo')
  return;
end

if ~isfield(data, 'cfg')
  % fieldtrip raw data structures are expected to have a cfg
  data.cfg = [];
end

hastrial   = isfield(data, 'trial');
hastime    = isfield(data, 'time');
hasfsample = isfield(data, 'fsample');

if ~hasfsample
  data.fsample = median(1./diff(data.time{1}));
end

if hastrial,
  ntrial = length(data.trial);
else
  ntrial = dimlength(data, 'rpt');
  if ~isfinite(ntrial) && strcmp(data.dimord(1:6), 'rpttap') && isfield(data, 'cumtapcnt'),
    ntrial = numel(data.cumtapcnt);
  elseif ~isfinite(ntrial)
    ntrial = 1;
  end
end

trl = ft_findcfg(data.cfg, 'trl');

nsmp = zeros(ntrial,1);
if hastrial,
  for i=1:ntrial
    nsmp(i) = size(data.trial{i}, 2);
  end
elseif ~isempty(trl)
  nsmp = trl(:,2) - trl(:,1) + 1;
end

if isempty(trl)
  warning('the data does not contain a trial definition, assuming that the trials are consecutive segments of a continuous recording');
  % construct a trial definition on the fly, assume that the trials are
  % consecutive segments of a continuous recording
  if ntrial==1,
    begsample = 1;
  else
    begsample = cat(1, 0, cumsum(nsmp(1:end-1))) + 1;
  end
  endsample = begsample + nsmp - 1;

  offset    = zeros(ntrial,1);
  if hastime,
    for i=1:ntrial
      offset(i) = time2offset(data.time{i}, data.fsample);
    end
  end
  trl = [begsample endsample offset];

elseif size(trl,1)~=ntrial
  warning('the trial definition in the configuration is inconsistent with the actual data');
  trl = [];
elseif nsmp~=(trl(:,2)-trl(:,1)+1)
  warning('the trial definition in the configuration is inconsistent with the actual data');
  trl = [];
end

if ~isfield(data, 'sampleinfo') && ~isempty(trl)
  data.sampleinfo = trl(:, 1:2);
elseif ~isfield(data, 'sampleinfo') && isempty(trl)
  warning('failed to create sampleinfo field');
end

if (~isfield(data, 'trialinfo') || isempty(data.trialinfo)) && ~isempty(trl) && size(trl, 2) > 3,
  data.trialinfo = trl(:, 4:end);
end

% if data is not raw then it does not make sense to keep the sampleinfo
if ~hastrial && isfield(data, 'sampleinfo')
  data = rmfield(data, 'sampleinfo');
end
