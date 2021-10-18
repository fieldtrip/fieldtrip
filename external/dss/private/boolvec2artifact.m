function artifact = boolvec2artifact(boolvec)

% BOOLVEC2ARTIFACT converts between two representations of events or trials.
%
% FieldTrip uses a number of representations for events that are conceptually very similar
%   event    = structure with type, value, sample, duration and offset
%   trl      = Nx3 numerical array with begsample, endsample, offset
%   trl      = table with 3 columns for begsample, endsample, offset
%   artifact = Nx2 numerical array with begsample, endsample
%   artifact = table with 2 columns for begsample, endsample
%   boolvec  = 1xNsamples boolean vector with a thresholded TTL/trigger sequence
%   boolvec  = MxNsamples boolean matrix with a thresholded TTL/trigger sequence
%
% If trl or artifact are represented as a MATLAB table, they can have additional
% columns. These additional columns have to be named and are not restricted to
% numerical values.
%
% This function makes an artifact definition from a Boolean vector, or a cell-array
% of artifact definitions from a Boolean matrix.
%
% See also ARTIFACT2BOOLVEC, ARTIFACT2EVENT, ARTIFACT2TRL, BOOLVEC2ARTIFACT, BOOLVEC2EVENT, BOOLVEC2TRL, EVENT2ARTIFACT, EVENT2BOOLVEC, EVENT2TRL, TRL2ARTIFACT, TRL2BOOLVEC, TRL2EVENT

% Copyright (C) 2009, Ingrid Nieuwenhuis
% Copyright (C) 2020, Robert Oostenveld
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

% represent it as a cell-array in case there are multiple channels
artifact = cell(1, size(boolvec,1));

for i=1:size(boolvec,1)
  tmp = diff([0 boolvec(i,:) 0]);
  begsample = find(tmp>0);
  endsample = find(tmp<0) - 1;
  artifact{i} = [begsample' endsample'];
end

if size(boolvec,1)==1
  % return it as a matrix, not as a cell-array of matrices
  artifact = artifact{1};
end
