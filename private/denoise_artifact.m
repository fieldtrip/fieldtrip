function [params, s_new] = denoise_artifact(params, s, state)

% DENOISE_ARTIFACT can be used for denoising source separation (DSS)
% during component analysis
%
% See also COMPONENTANALYSIS

% Copyright (C) 2005, Robert Oostenveld
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

if nargin<3 || ~isstruct(state)
  params.name = 'denoise_artifact';
  return;
end

if ~isfield(params, 'concat')
  % the artifact matrix should be relative to the concatenated data
  % and not to the samples in the data file
  trl      = params.trl;
  artifact = params.artifact;
  concat   = [];
  trlnumsmp = trl(:,2) - trl(:,1) + 1;   % count the number of samples in each trial
  trlnumsmp = [0; cumsum(trlnumsmp)];    % compute the cumulative number of samples
  for i=1:size(artifact,1)
    % find the trial in which this artifact lies
    trlnum = find(artifact(i,1)>=trl(:,1) & artifact(i,2)<=trl(:,2));
    if length(trlnum)==0
      fprintf('artifact %d is not completely in one trial\n', i);
    elseif length(trlnum)==1
      fprintf('artifact %d is completely in trial %d\n', i, trlnum);
      artbeg = artifact(i,1);
      artend = artifact(i,2);
      trlbeg = trl(trlnum,1);
      trlend = trl(trlnum,2);
      dum1 = artbeg - trlbeg + trlnumsmp(trlnum) + 1;
      dum2 = artend - trlbeg + trlnumsmp(trlnum) + 1;
      concat(end+1,:) = [dum1 dum2];
    else
      fprintf('artifact %d seems to overlap with multiple trials\n', i);
    end
  end
  fprintf('%d trial segments\n', size(trl,1));
  fprintf('%d artifact segments in input\n', size(artifact,1));
  fprintf('%d artifact segments in output\n', size(concat,1));
  % remember the artifact sample numbers in the concatinated data
  params.concat = concat;
end

% the artifacts are expressed in samples w.r.t. the concatenated data
artifact = params.concat;

begsmp = artifact(:,1);
endsmp = artifact(:,2);

% FIXME, this assumes that all artifacts are equally long
numsmp = endsmp(1) - begsmp(1) + 1;
numsrc = size(s,1);
numart = size(artifact,1);
sum    = zeros(numsrc,numsmp);

for i=1:numart
  % cut data segments around the ECG trigger
  sum = sum + s(:,begsmp(i):endsmp(i));
end

% averaging is not strictly neccessary, but conceptually nicer
sum = sum ./ numart;

s_new = zeros(size(s));
for i=1:numart
  % put the average data segments back around the ECG trigger
  s_new(:,begsmp(i):endsmp(i)) = sum;
end

