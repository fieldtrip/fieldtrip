function data = fixsampleinfo(data)

% FIXSAMPLEINFO checks for the existence of a sampleinfo field in the
% provided data structure. If present, nothing is done; if absent,
% fixsampleinfo attempts to reconstruct a sampleinfo based on either an
% trl-matrix present in the cfg-tree, or by just assuming the trials are
% segments of a continuous recording.

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

% check whether we're dealing with a timelock structure that has trials
istimelock = hastime && hastrial && ~iscell(data.trial) && ~iscell(data.time);

if ~hasfsample && hastime
  if istimelock
    data.fsample = median(1./diff(data.time));
  else
    data.fsample = median(1./diff(data.time{1}));
  end
end

if hastrial
  if istimelock
    ntrial = size(data.trial,1);
  else
    ntrial = numel(data.trial);
  end
else
  ntrial = dimlength(data, 'rpt');
  if ~isfinite(ntrial) && strcmp(data.dimord(1:6), 'rpttap') && isfield(data, 'cumtapcnt'),
    ntrial = numel(data.cumtapcnt);
  elseif ~isfinite(ntrial)
    ntrial = 1;
  end
end

trl = ft_findcfg(data.cfg, 'trl');

if istimelock
  nsmp = ones(ntrial,1) .* size(data.trial,3);
else
  nsmp = zeros(ntrial,1);
  if hastrial
    for i=1:ntrial
      nsmp(i) = size(data.trial{i}, 2);
    end
  elseif ~isempty(trl)
    nsmp = trl(:,2) - trl(:,1) + 1;
  end
end

if isempty(trl)
  warning_once('the data does not contain a trial definition');
elseif ~isempty(trl) && size(trl,1)~=numel(nsmp)
  warning_once('the trial definition in the configuration is inconsistent with the actual data');
  trl = [];
elseif size(trl,1)~=ntrial
  warning_once('the trial definition in the configuration is inconsistent with the actual data');
  trl = [];
elseif nsmp~=(trl(:,2)-trl(:,1)+1)
  warning_once('the trial definition in the configuration is inconsistent with the actual data');
  trl = [];
end

if isempty(trl) || ~all(nsmp==trl(:,2)-trl(:,1)+1)
  warning_once('reconstructing sampleinfo by assuming that the trials are consecutive segments of a continuous recording');
  % construct a trial definition on the fly, assume that the trials are
  % consecutive segments of a continuous recording
  if ntrial==1,
    begsample = 1;
  else
    begsample = cat(1, 0, cumsum(nsmp(1:end-1))) + 1;
  end
  endsample = begsample + nsmp - 1;

  if istimelock
    offset = ones(ntrial,1) .* time2offset(data.time, data.fsample);
  else
    offset    = zeros(ntrial,1);
    if hastime,
      for i=1:ntrial
        offset(i) = time2offset(data.time{i}, data.fsample);
      end
    end
  end
  trl = [begsample endsample offset];
end

if ~isfield(data, 'sampleinfo') && ~isempty(trl)
  data.sampleinfo = trl(:, 1:2);
elseif ~isfield(data, 'sampleinfo') && isempty(trl)
  % this is probably an unreachable statement
  warning_once('failed to create sampleinfo field');
end

if (~isfield(data, 'trialinfo') || isempty(data.trialinfo)) && ~isempty(trl) && size(trl, 2) > 3,
  data.trialinfo = trl(:, 4:end);
end

% if data does not have repetitions (i.e. trials) then it does not make sense to keep the sampleinfo
if ~hastrial && isfield(data, 'sampleinfo')
  data = rmfield(data, 'sampleinfo');
end
