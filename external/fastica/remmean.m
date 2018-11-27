function [newVectors, meanValue] = remmean(vectors)
%REMMEAN - remove the mean from vectors
%
% [newVectors, meanValue] = remmean(vectors);
%
% Removes the mean of row vectors.
% Returns the new vectors and the mean.
%
% This function is needed by FASTICA and FASTICAG

% @(#)$Id$

%newVectors = zeros (size (vectors));
meanValue = mean (vectors,2);
if ~iscell(vectors)
  newVectors = vectors - meanValue * ones (1,size (vectors, 2));
else
  newVectors = cellvecadd(vectors,-meanValue);
end