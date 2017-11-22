function [numA, numB, indA, indB] = spikesort(numA, numB, varargin)

% SPIKESORT uses a variation on the cocktail sort algorithm in combination
% with a city block distance to achieve N-D trial pairing between spike
% counts. The sorting is not guaranteed to result in the optimal pairing. A
% linear pre-sorting algorithm is used to create good initial starting
% positions.
%
% The goal of this function is to achieve optimal trial-pairing prior to
% stratifying the spike numbers in two datasets by random removal of some
% spikes in the trial and channel with the largest numnber of spikes.
% Pre-sorting based on the city-block distance between the spike count
% ensures that as few spikes as possible are lost.
%
% Use as
%   [srtA, srtB, indA, indB] = spikesort(numA, numB, ...)
%
% Optional arguments should be specified as key-value pairs and can include
%   'presort'  number representing the column, 'rowwise' or 'global'
%
% Example
%   numA = reshape(randperm(100*3), 100, 3);
%   numB = reshape(randperm(100*3), 100, 3);
%   [srtA, srtB, indA, indB] = spikesort(numA, numB);
%   % check that the order is correct, the following should be zero
%   numA(indA,:) - srtA
%   numB(indB,:) - srtB
%
% See also COCKTAILSORT

% Copyright (C) 2007, Robert Oostenveld
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

% this can be used for printing detailled user feedback
fb = false;

% get the options
presort = ft_getopt(varargin, 'presort');

if any(size(numA)~=size(numB))
  ft_error('input dimensions should be the same');
end

bottom      = 1;
top         = size(numA,1);
swapped     = true;
dist        = zeros(1,top);  % this will hold the distance between the trials in each pair
distswap    = zeros(1,2);    % this will hold the distance for two trials after swapping

% this is to keep track of the row-ordering during sorting
sel = 1:size(numA,2);
numA(:,end+1) = 1:top;
numB(:,end+1) = 1:top;

% compute the initial distance between the trial pairs using city block distance metric
for i=1:top
  dist(i) = sum(abs(numA(i,sel)-numB(i,sel)));
end

if fb
  fprintf('initial cost        = %d\n', sum(dist));
end

if isnumeric(presort)
  % start by pre-sorting on the first column only
  [dum, indA] = sort(numA(:,presort));
  [dum, indB] = sort(numB(:,presort));
  numA = numA(indA,:);
  numB = numB(indB,:);
elseif strcmp(presort, 'rowwise')
  full = cityblock(numA(:,sel), numB(:,sel));
  link = zeros(1,top);
  for i=1:top
    d = full(i,:);
    d(link(1:(i-1))) = inf;
    [m, ind] = min(d);
    link(i) = ind;
  end
  indA = 1:top;
  indB = link;
  numB = numB(link,:);
elseif strcmp(presort, 'global')
  full = cityblock(numA(:,sel), numB(:,sel));
  link = zeros(1,top);
  while ~all(link)
    [m, i] = min(full(:));
    [j, k] = ind2sub(size(full), i);
    full(j,:) = inf;
    full(:,k) = inf;
    link(j) = k;
  end
  indA = 1:top;
  indB = link;
  numB = numB(link,:);
end

% compute the initial distance between the trial pairs
% using city block distance metric
for i=1:top
  dist(i) = sum(abs(numA(i,sel)-numB(i,sel)));
end

if fb
  fprintf('cost after pre-sort = %d\n', sum(dist));
end

while swapped
  swapped = false;

  for i=bottom:(top-1)

    % compute the distances between the trials after swapping
    % using city block distance metric
    distswap(1) = sum(abs(numA(i,sel)-numB(i+1,sel)));
    distswap(2) = sum(abs(numA(i+1,sel)-numB(i,sel)));

    costNow = sum(dist([i i+1]));
    costSwp = sum(distswap);

    if costNow>costSwp                        % test whether the two elements are in the correct order
      numB([i i+1], :) = numB([i+1 i], :);    % let the two elements change places
      dist([i i+1]) = distswap;               % update the distance vector
      swapped = true;
    end
  end

  % decreases `top` because the element with the largest value in the unsorted
  % part of the list is now on the position top
  top = top - 1;

  for i=top:-1:(bottom+1)

    % compute the distances between the trials after swapping
    % using city block distance metric
    distswap(1) = sum(abs(numA(i,sel)-numB(i-1,sel)));
    distswap(2) = sum(abs(numA(i-1,sel)-numB(i,sel)));
    costNow = sum(dist([i i-1]));
    costSwp = sum(distswap);

    if costNow>costSwp
      numB([i i-1], :) = numB([i-1 i], :);    % let the two elements change places
      dist([i i-1]) = distswap;               % update the distance vector
      swapped = true;
    end
  end

  % increases `bottom` because the element with the smallest value in the unsorted
  % part of the list is now on the position bottom
  bottom = bottom + 1;

end % while swapped

if fb
  fprintf('final cost          = %d\n', sum(dist));
end

indA = numA(:,end);
numA = numA(:,sel);
indB = numB(:,end);
numB = numB(:,sel);
end % function spikesort

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = cityblock(a, b)
d = zeros(size(a,1), size(b,1));
for i=1:size(a,1)
  for j=1:size(b,1)
    d(i,j) = sum(abs(a(i,:)-b(j,:)));
  end
end
end % function cityblock
