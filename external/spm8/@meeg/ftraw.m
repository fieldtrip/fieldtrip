function raw = ftraw(this, chanind, timeind, trialind)
% Method for converting meeg object to Fieldtrip raw struct
% FORMAT raw  =  ftraw(this, chanind, timeind, trialind)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: ftraw.m 5063 2012-11-16 12:10:41Z vladimir $

if ~isequal(transformtype(this), 'time')
    raw = fttimelock(this, chanind, timeind, trialind);
    return;
end

% chanind == 0 is accepted for backward compatibility
if nargin < 2 || isempty(chanind) || chanind == 0
    chanind = 1:nchannels(this);
end

if nargin < 3 || isempty(timeind)
    timeind = 1:nsamples(this);
end

if nargin < 4 || isempty(trialind)
    trialind = 1:ntrials(this);
end

raw = [];

raw.label   = chanlabels(this, chanind)';

raw.trial   = cell(1, length(trialind));

for i =  1:length(trialind) 
    raw.trial{i} = this.data.y(chanind, timeind, trialind(i));
end

raw.time = repmat({time(this, timeind)}, 1, length(trialind));

clist      =  condlist(this);

condlabels = conditions(this, trialind);

raw.trialinfo = 0*trialind;
for k = 1:numel(clist)
  fprintf('mapping condition label "%s" to condition code %d\n', clist{k}, k);
  sel = strcmp(clist{k}, condlabels);
  raw.trialinfo(sel) = k;
end

if ~isempty(sensors(this, 'MEG'))
    raw.grad = sensors(this, 'MEG');
end

if ~isempty(sensors(this, 'EEG'))
    raw.elec = sensors(this, 'EEG');
end

if isfield(this.other, 'origheader')
    raw.hdr = this.other.origheader;
end

onsets = trialonset(this, trialind);

if all(onsets>0)
    onsets = round(onsets(:)*fsample(this));
    raw.sampleinfo = [onsets+timeind(1) onsets+timeind(end)]-1;
end

