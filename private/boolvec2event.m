function event = boolvec2event(boolvec, varargin)

% BOOLVEC2EVENT converts between two representations of events or trials.
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

% get the optional arguments or set defaults
type  = ft_getopt(varargin, 'type', repmat({'trigger'}, 1, size(boolvec, 1)));
value = ft_getopt(varargin, 'value', {}); % by default it will take the value from the vector or matrix

if ~iscell(type)
  type = {type};
end

if ~iscell(value) 
  value = {value};
end

event = [];

for i=1:size(boolvec,1)
  tmp = diff([0 boolvec(i,:) 0]);
  begsample = find(tmp>0);
  endsample = find(tmp<0) - 1;
  for j=1:numel(begsample)
    event(end+1).type     = type{i};
    if ~isempty(value)
      event(end  ).value  = value{i};
    else
      event(end  ).value  = boolvec(i, begsample(j));
    end
    event(end  ).sample   = begsample(j);
    event(end  ).duration = endsample(j) - begsample(j) + 1;
    event(end  ).offset   = 0;
  end
end
