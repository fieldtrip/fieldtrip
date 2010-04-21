function data = fixtrialdef(data)

% FIXTRIALDEF adds a trl matrix to the raw data configuration
% which is constructed on the fly, assuming that the trials are
% consecutive segments of a continuous recording

% Copyright (C) 2009, Robert Oostenveld

if ~isfield(data, 'cfg')
  % fieldtrip raw data structures are expected to have a cfg
  data.cfg = [];
end

ntrial = length(data.trial);
trl    = findcfg(data.cfg, 'trl');

nsmp = zeros(ntrial,1);
for i=1:ntrial
  nsmp(i) = size(data.trial{i}, 2);
end

if isempty(trl)
  warning('the data does not contain a trial definition, assuming that the trials are consecutive segments of a continuous recording');
  % construct a trial definition on the fly, assume that the trials are
  % consecutive segments of a continuous recording
  begsample = cat(1, 0, cumsum(nsmp(1:end-1))) + 1;
  endsample = begsample + nsmp - 1;
  offset    = zeros(ntrial,1);
  for i=1:ntrial
    offset(i) = time2offset(data.time{i}, data.fsample);
  end
  data.cfg.trl = [begsample endsample offset];

elseif size(trl,1)~=ntrial
  error('the trial definition in the configuration is inconsistent with the actual data');
elseif nsmp~=(trl(:,2)-trl(:,1)+1)
  error('the trial definition in the configuration is inconsistent with the actual data');
end
