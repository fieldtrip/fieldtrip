function [bp, x] = binomialprob(pobs, alpha, subjratio)

% BINOMIALPROB computes the probability of observing a significant effect
% in multiple tests. It allows you to test questions like "How likely
% is it that there is a significant effect at this time-frequency point
% for 8 out of 10 subjects, given that the probability of observing a
% significant effect in a given subject is 5%"
% 
% Use as
%    [bprob] = binomialprob(prob, alpha)
% where
%   prob   is a Nvoxel X Nsubject matrix with the single-subject probability
%   alpha  is the probability of observing a significant voxel
%
% The function also has more advanced functionality, please read the code 
% if you are interested.
%
% See also BINOPDF, BINOCDF

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

% determine the number of subjects
[M, N] = size(pobs);

if nargin<2
  alpha = [];
end

if nargin<3
  subjratio = [];
end

if ~isempty(subjratio)
  % threshold the statistical maps per subject to obtain the desired ratio
  for i=1:N
    s = sort(pobs(:,i));           % sort the values over voxels
    t = round((1-subjratio) * M);  % determine the index of the threshold
    a = (pobs(:,i)>=s(t));         % determine the voxels that exceed the threshold
    pobs(:,i) = a;                 % assign the thresholded statistic
  end
end

% determine whether the single subject statistical maps have already been thresholded
isthresh = all(pobs(:)==0 | pobs(:)==1);

% determine whether the probability for the binomial distribution is specified
isalpha  = ~isempty(alpha);

if      isthresh &&  isalpha
  % the probability of observing a significant voxel is specified
  p = alpha;
  % count the number of subjects in which each voxel is significant
  x = sum(pobs, 2);
elseif  isthresh && ~isalpha
  % estimate the probability from the ratio of thresholded voxels
  p = sum(pobs(:))/length(pobs(:));
  % count the number of subjects in which each voxel is significant
  x = sum(pobs, 2);
elseif ~isthresh &&  isalpha
  % the probability of observing a significant voxel is specified
  p = alpha;
  % threshold the single subject probability maps at the alpha level
  x = sum(pobs<=alpha, 2);
elseif ~isthresh && ~isalpha
  error('can only determine alpha automatically from thresholded statistical maps');
end

% this uses MATLAB stats toolbox
bp = 1 - binocdf(x, N, p);

