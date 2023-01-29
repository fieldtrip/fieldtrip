function data = fixsampleinfo(data)

% FIXSAMPLEINFO checks for the existence of a sampleinfo and trialinfo field in the
% provided raw or timelock data structure. If present, nothing is done; if absent,
% this function attempts to reconstruct them based on either an trl-matrix present in
% the cfg-tree, or by just assuming the trials are segments of a continuous
% recording.
%
% See also FT_DATATYPE_RAW, FT_DATATYPE_TIMELOCK

% Copyright (C) 2009-2020, Robert Oostenveld and Jan-Mathijs Schoffelen

if isfield(data, 'sampleinfo')
  % sampleinfo is present (as requested), trialinfo is optional
  return
end


hastrial   = isfield(data, 'trial');
hastime    = isfield(data, 'time');
hasfsample = isfield(data, 'fsample');

% check whether we're dealing with a timelock structure that has trials
istimelock = hastime && hastrial && ~iscell(data.trial) && ~iscell(data.time);
israw      = hastime && hastrial &&  iscell(data.trial) &&  iscell(data.time);

% if the data does not have repetitions (i.e. trials or subjects), then it does not make sense to have these fields
if ~hastrial && (isfield(data, 'sampleinfo') || isfield(data, 'trialinfo'))
  data = removefields(data, {'sampleinfo', 'trialinfo'});
  return
end

if ~hasfsample && hastime
  if israw
    data.fsample = median(1./diff(data.time{1}));
  elseif hastime
    data.fsample = median(1./diff(data.time));
  end
end

if hastrial
  if israw
    ntrial = numel(data.trial);
  elseif istimelock
    dimord = getdimord(data, 'trial');
    switch dimord
      case {'rpt_chan_time', 'subj_chan_time'}
        ntrial = size(data.trial, 1);
      case {'chan_time'}
        ntrial = 1;
      otherwise
        ft_error('unexpected dimord')
    end % switcch
  end
else
  ntrial = dimlength(data, 'rpt');
  if ~isfinite(ntrial) && startsWith(data.dimord, 'rpttap') && isfield(data, 'cumtapcnt')
    ntrial = numel(data.cumtapcnt);
  elseif ~isfinite(ntrial)
    ntrial = 1;
  end
end

if isfield(data, 'cfg')
  % this might find and return something, but could also return empty
  trl = ft_findcfg(data.cfg, 'trl');
else
  trl = [];
end

if istable(trl)
  % the subsequent code requires this to be an array
  trl = table2array(trl(:,1:3));
end

if israw
  nsmp = zeros(ntrial,1);
  if israw
    for i=1:ntrial
      nsmp(i) = size(data.trial{i}, 2);
    end
  elseif ~isempty(trl)
    nsmp = trl(:,2) - trl(:,1) + 1;
  end
elseif istimelock
  nsmp = ones(ntrial,1) .* size(data.trial, ndims(data.trial));
elseif hastime
  nsmp = ones(ntrial,1) .* length(data.time);
end

if isempty(trl)
  ft_warning('the data does not contain sampleinfo');
elseif ~isempty(trl) && size(trl,1)~=numel(nsmp)
  ft_warning('sampleinfo in the configuration is inconsistent with the actual data');
  trl = [];
elseif size(trl,1)~=ntrial
  ft_warning('sampleinfo in the configuration is inconsistent with the actual data');
  trl = [];
elseif nsmp~=(trl(:,2)-trl(:,1)+1)
  ft_warning('sampleinfo in the configuration is inconsistent with the actual data');
  trl = [];
end

if isempty(trl) || ~all(nsmp==trl(:,2)-trl(:,1)+1)
  ft_warning('reconstructing sampleinfo by assuming that the trials are consecutive segments of a continuous recording');
  % construct the trial definition on the fly
  if ntrial==1
    begsample = 1;
  else
    begsample = cat(1, 0, cumsum(nsmp(1:end-1))) + 1;
  end
  endsample = begsample + nsmp - 1;
  
  if israw
    offset = zeros(ntrial,1);
    if hastime
      for i=1:ntrial
        offset(i) = time2offset(data.time{i}, data.fsample);
      end
    end
  elseif hastime
    offset = ones(ntrial,1) .* time2offset(data.time, data.fsample);
  end
  trl = [begsample endsample offset];
end

if ~isfield(data, 'sampleinfo') && ~isempty(trl)
  data.sampleinfo = trl(:, 1:2);
elseif ~isfield(data, 'sampleinfo') && isempty(trl)
  % this is probably an unreachable statement
  ft_warning('failed to create sampleinfo');
end

if (~isfield(data, 'trialinfo') || isempty(data.trialinfo)) && size(trl, 2) > 3
  % update the trialinfo using the additional columns from the trial definition
  data.trialinfo = trl(:, 4:end);
end
