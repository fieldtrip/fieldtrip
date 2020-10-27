function trl = sampleinfo2trl(data)

% SAMPLEINFO2TRL constructs the trial definition from the sampleinfo, the time axes
% and optionally from the trialinfo
%
% Use as
%   trl = sampleinfo2trl(data)
%
% See also ARTIFACT2BOOLVEC, ARTIFACT2EVENT, ARTIFACT2TRL, BOOLVEC2ARTIFACT, BOOLVEC2EVENT, BOOLVEC2TRL, EVENT2ARTIFACT, EVENT2BOOLVEC, EVENT2TRL, TRL2ARTIFACT, TRL2BOOLVEC, TRL2EVENT

% get the begin and end sample of each trial
begsample = data.sampleinfo(:,1);
endsample = data.sampleinfo(:,2);

% recreate the offset
offset = zeros(numel(data.trial), 1);
for i=1:numel(data.trial)
  offset(i) = time2offset(data.time{i}, data.fsample);
end

if isfield(data, 'trialinfo') && istable(data.trialinfo)
  trl = table(begsample, endsample, offset);
  trl = horzcat(trl, data.trialinfo);
elseif isfield(data, 'trialinfo') && isnumeric(data.trialinfo)
  trl = [begsample endsample offset data.trialinfo];
else
  trl = [begsample endsample offset];
end
