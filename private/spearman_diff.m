function s = spearman_diff(cfg, dat, design)

%SPEARMAN_DIFF computes the difference in rank-correlation 
%coefficient between two variables, between conditions '1'
%and '2' as they are labeled by design 

%The input-data should be formatted as follows:
% first dimension : signals (or signal-combinations)
% second dimension: repetitions
% third dimension : frequencies (optional), or the two signals which will be correlated
% fourth dimension (optional): the two signals which will be correlated
% the last dimension should have length two, since this dimension contains the two variables 
% that are to be rank-correlated

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

if ~isfield(cfg, 'factor') & prod(size(design)) ~= max(size(design)),
  error('cannot determine the labeling of the trials');
elseif ~isfield(cfg, 'factor')
  cfg.ivar = 1;
end

nsgn = size(dat,1);
nrpt = size(dat,2);
if length(size(dat))==3,
  nfrq = 1;
  n    = size(dat,3);
else
  nfrq = size(dat,3);
  n    = size(dat,4);
end

if n ~= 2,
  error('the last dimension of the input should be 2');
end

selA = find(design(cfg.ivar, :) == 1);
selB = find(design(cfg.ivar, :) == 2);

%if length(selA) ~= length(selB)
%  error('inappropriate design');
%end
%NONSENSE

for k = 1:nsgn
  for j = 1:nfrq
     datA = squeeze(dat(k, selA, j, :)); %datA = trials x 2
     datB = squeeze(dat(k, selB, j, :)); %datB = trials x 2
     
     [srtA, indA]     = sort(datA); 
     datA(indA(:,1),1) = [1:size(datA,1)]';
     datA(indA(:,2),2) = [1:size(datA,1)]';
     
     [srtB, indB]     = sort(datB); 
     datB(indB(:,1),1) = [1:size(datB,1)]';
     datB(indB(:,2),2) = [1:size(datB,1)]';

     denomA     = size(datA,1) * (size(datA,1)^2-1) / 6;
     rccA(k, j) = 1 - sum(diff(datA, [], 2).^2) / denomA;

     denomB     = size(datB,1) * (size(datB,1)^2-1) / 6;
     rccB(k, j) = 1 - sum(diff(datB, [], 2).^2) / denomB;
  end
end

s = rccA - rccB;

