function trl = event2trl(event)

% EVENT2TRL converts between two representations of events or trials.
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

type     = {event.type}';     % these are strings
value    = {event.value}';    % can be string or number
sample   = [event.sample]';   % these are numbers
offset   = {event.offset}';   % can be empty or a number
duration = {event.duration}'; % can be empty or a number

% do some sanity checks
assert(~any(cellfun(@isempty, type)));
assert(~any(cellfun(@isempty, value)));
assert(length(sample)==numel(event));

% set default to 0 and convert to vector
offset(cellfun(@isempty, offset)) = {0};
offset = [offset{:}]';

% set default to 1 and convert to vector
duration(cellfun(@isempty, duration)) = {1};
duration = [duration{:}]';

% set the minimum duration to one sample
duration(duration<1) = 1;

if all(cellfun(@isnumeric, value))
  value = [value{:}]';
end

begsample = sample;
endsample = begsample + duration - 1;

if length(unique(type))==1
  if iscell(value) && all(ischar(value)) && length(unique(value))==1
    % don't add the type or value
    trl = [begsample endsample offset];
  elseif all(isnumeric(value))
    % add the value as a number
    trl = [begsample endsample offset value];
  else
    % add the value as a string, this requires the trl to be a MATLAB table
    trl = table(begsample, endsample, offset);
    trl.value = value;
  end
else
  if iscell(value) && all(ischar(value)) && length(unique(value))==1
    % ignore the value but add the type as a string, this requires the trl to be a MATLAB table
    trl = table(begsample, endsample, offset);
    trl.type = type;
  else
    % add both the type and value as a string, this requires the trl to be a MATLAB table
    trl = table(begsample, endsample, offset);
    trl.type = type;
    trl.value = value;
  end
end
