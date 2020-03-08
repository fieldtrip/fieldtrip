function trl = sampleinfo2trl(data)

% SAMPLEINFO2TRL constructs the trial definition from the sampleinfo, the time axes
% and optionally from the trialinfo

% get the begin and end sample of each trial
trl = zeros(numel(data.trial), 3);
trl(:,[1 2]) = data.sampleinfo;

% recreate the offset
for ntrl = 1:numel(data.trial)
  trl(ntrl,3) = time2offset(data.time{ntrl}, data.fsample);
end

if isfield(data, 'trialinfo')
  if istable(data.trialinfo)
    % convert table into normal array
    trl = [trl table2array(data.trialinfo)];
  else
    trl = [trl data.trialinfo];
  end
end
