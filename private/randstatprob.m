function p = randstatprob(randobs, realobs, tail, correctm)

% RANDSTATPROB computes the non-parametric probability of the observed
% value under the assumption that the random observations are equally
% probable under the null hypothesis.
%
% Use as
%   p = randstatprob(randobs, realobs, tail, correctm)
% where
%   randobs  = Nvox x Nrnd
%   realobs  = Nvox x 1, or Nvox x Nobs (for multiple observations)
%   tail     =  0 for two-sided test
%   tail     =  1 for one-sided test with realobs>=randobs
%   tail     = -1 for one-sided test with realobs<=randobs
%   correctm =  0 do not correct for multiple comparisons
%               1 correct for multiple comparisons using the maximum statistic
%               2 correct for multiple comparisons using ordered statistics
%
% Each row of the input data contains all the (real or randomized)
% observations in one voxel. Multiple comparison can be performed by
% creating a reference distribution based on the minimum or maximum
% of all voxels for each randomization.

% Copyright (C) 2004-2005, Robert Oostenveld
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

if nargin<3
  tail = 0;
end

if nargin<4
  correctm = 0;
end

Nvox = size(randobs,1);
Nrnd = size(randobs,2);
Nobs = size(realobs,2);
if size(realobs,1)~=Nvox
  error('dimensions of input arguments does not match');
end

p = zeros(size(realobs));

if correctm==1
  % the distribution of the min/max statistic is the same for each voxel
  mindist = min(randobs, [], 1);
  maxdist = max(randobs, [], 1);
elseif correctm==2
  % sort the observed and the randomized statistics for the ordered statistics
  indx = zeros(size(realobs));
  for i=1:Nobs
    [realobs(:,i), indx(:,i)] = sort(realobs(:,i));
  end
  for i=1:Nrnd
    randobs(:,i) = sort(randobs(:,i));
  end
end

for i=1:Nobs
  if correctm==0 || correctm==2
    % do not apply multiple comparison correction or use ordered statistics
    switch tail
      case 0
        plo = sum(randobs <= repmat(realobs(:,i),1,Nrnd), 2)./Nrnd;
        phi = sum(randobs >= repmat(realobs(:,i),1,Nrnd), 2)./Nrnd;
        p(:,i) = min(plo,phi);

      case 1
        p(:,i) = sum(randobs >= repmat(realobs(:,i),1,Nrnd), 2)./Nrnd;

      case -1
        p(:,i) = sum(randobs <= repmat(realobs(:,i),1,Nrnd), 2)./Nrnd;

      otherwise
        error('incorrect specification of tail');
    end
  elseif correctm==1
    % apply multiple comparison correction using the maximum statistic
    switch tail
      case 0
        % the observation can be either in the tail with the smallest or the tail with the largest values
        plo = sum(repmat(mindist,Nvox,1) <= repmat(realobs(:,i),1,Nrnd), 2)./Nrnd;
        phi = sum(repmat(maxdist,Nvox,1) >= repmat(realobs(:,i),1,Nrnd), 2)./Nrnd;
        p(:,i) = min(plo,phi);

      case 1
        p(:,i) = sum(repmat(maxdist,Nvox,1) >= repmat(realobs(:,i),1,Nrnd), 2)./Nrnd;

      case -1
        p(:,i) = sum(repmat(mindist,Nvox,1) <= repmat(realobs(:,i),1,Nrnd), 2)./Nrnd;

      otherwise
        error('incorrect specification of tail');
    end
  end
end

if correctm==2
  % undo the sorting that was needed for the ordered statistics
  for i=1:Nobs
    [dum, unsort] = sort(indx(:,i));
    p(:,i) = p(unsort,i);
  end
end

