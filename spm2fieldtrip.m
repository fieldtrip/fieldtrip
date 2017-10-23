function data = spm2fieldtrip(D)

% SPM2FIELDTRIP converts an SPM8 meeg object into a FieldTrip raw data structure
%
% Use as
%   data = spm2fieldtrip(D)
% where D is the SPM meeg object which you can load in with SPM_EEG_LOAD
% and where data is a FieldTrip raw data structure as if it were returned
% by FT_PREPROCESSING.
%
% See also FT_PREPROCESSING, SPM_EEG_LOAD

if ~ft_hastoolbox('spm8up')
  % it should be version spm8 or higher, since spm99, spm2 and spm5 did not yet the "meeg" object
  ft_error('this requires SPM8 or later to be on your MATLAB path');
end

if ~isa(D, 'meeg')
  ft_error('this requires an SPM "meeg" object as input')
end

% this is how SPM8 represents it
data = D.ftraw(0);

clist      = D.condlist;
conditions = D.conditions;

data.trialinfo = zeros(D.ntrials,1);
for k = 1:numel(clist)
  fprintf('mapping condition label "%s" to condition code %d\n', clist{k}, k);
  sel=strcmp(clist{k}, conditions);
  data.trialinfo(sel) = k;
end

% FIXME the following is not correct
%
% data.sampleinfo = zeros(ntrials,2);
% for i=1:ntrials
%   data.sampleinfo(i,1) = D.indsample(i);
%   data.sampleinfo(i,2) = D.indsample(i) + D.nsamples;
% end
