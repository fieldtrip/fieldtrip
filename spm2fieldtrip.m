function [data] = spm2fieldtrip(D)

% SPM2FIELDTRIP converts an SPM8 meeg object into a FieldTrip raw data structure
%
% Use as
%   data = spm2fieldtrip(D)
% where D is the SPM meeg object which you can load in with SPM_EEG_LOAD
% and where data is a FieldTrip raw data structure as if it were returned
% by FT_PREPROCESSING.
%
% See also FT_PREPROCESSING, SPM_EEG_LOAD

% Copyright (C) 2011, Donders Centre for Cognitive Neuroimaging
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

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
