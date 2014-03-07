function [data] = spm2fieldtrip(D)

% SPM2FIELDTRIP converts an SPM8 meeg object into a FieldTrip raw data structure
%
% Use as
%   data = spm2fieldtrip(D)
% where D is the SPM8 meeg object which you can load in with SPM_EEG_LOAD
% and where data is a FieldTrip raw data structure as if it were returned
% by FT_PREPROCESSING.
%
% See also FT_PREPROCESSING, SPM_EEG_LOAD

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble provenance

if ~ft_hastoolbox('spm8')
  error('this requires a full version of SPM8 on your MATLAB path');
end

if ~isa(D, 'meeg')
  error('this requires an SPM8 "meeg" object as input')
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

% create empty cfg-structure in order for the ft_postamble to work. It is
% of no further consequence
cfg = [];

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble history data

