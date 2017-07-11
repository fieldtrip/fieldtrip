function event = ft_filter_event(event, varargin)

% FT_FILTER_EVENT does what its name implies
%
% Use as
%   event = ft_filter_event(event, ...)
%
% The optional arguments should come in key-value pairs and determine the
% filter characteristics:
%   type         = cell-array with strings
%   value        = numeric array
%   sample       = numeric array
%   timestamp    = numeric array
%   offset       = numeric array
%   duration     = numeric array
%   minsample    = value
%   maxsample    = value
%   minduration  = value
%   maxduration  = value
%   mintimestamp = value
%   maxtimestamp = value
%   minnumber    = value, applies only if event.number is present
%   maxnmumber   = value, applies only if event.number is present
%
% See also FT_READ_EVENT, FT_WRITE_EVENT

% Copyright (C) 2007-2010 Robert Oostenveld
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

% get the optional input arguments
type         = ft_getopt(varargin, 'type');
value        = ft_getopt(varargin, 'value');
sample       = ft_getopt(varargin, 'sample');
timestamp    = ft_getopt(varargin, 'timestamp');
offset       = ft_getopt(varargin, 'offset');
duration     = ft_getopt(varargin, 'duration');

% the numeric fields can also be filtered on a range
minsample    = ft_getopt(varargin, 'minsample');
maxsample    = ft_getopt(varargin, 'maxsample');
minduration  = ft_getopt(varargin, 'minduration');
maxduration  = ft_getopt(varargin, 'maxduration');
mintimestamp = ft_getopt(varargin, 'mintimestamp');
maxtimestamp = ft_getopt(varargin, 'maxtimestamp');
minnumber    = ft_getopt(varargin, 'minnumber');
maxnumber    = ft_getopt(varargin, 'maxnumber');

if ~isempty(type) && ~iscell(type)
  % this can be specified as string or as cell-array, convert to cell-array
  type = {type};
end

% determine which filters to apply
testtype         = ~isempty(type)         && isfield(event, 'type');
testvalue        = ~isempty(value)        && isfield(event, 'value');
testsample       = ~isempty(sample)       && isfield(event, 'sample');
testtimestamp    = ~isempty(timestamp)    && isfield(event, 'timestamp');
testoffset       = ~isempty(offset)       && isfield(event, 'offset');
testduration     = ~isempty(duration)     && isfield(event, 'duration');
testminsample    = ~isempty(minsample)    && isfield(event, 'sample');
testmaxsample    = ~isempty(maxsample)    && isfield(event, 'sample');
testminduration  = ~isempty(minduration)  && isfield(event, 'duration');
testmaxduration  = ~isempty(maxduration)  && isfield(event, 'duration');
testmintimestamp = ~isempty(mintimestamp) && isfield(event, 'timestamp');
testmaxtimestamp = ~isempty(maxtimestamp) && isfield(event, 'timestamp');
testminnumber    = ~isempty(minnumber)    && isfield(event, 'number');
testmaxnumber    = ~isempty(maxnumber)    && isfield(event, 'number');

if (~isempty(minnumber) || ~isempty(maxnumber)) && ~isfield(event, 'number')
  ft_warning('the events are not numbered, assuming that the order corresponds to the original stream sequence');
  for i=1:length(event)
    event(i).number = i;
  end
  testminnumber    = ~isempty(minnumber);
  testmaxnumber    = ~isempty(maxnumber);
end

% apply the filters
sel = true(length(event),1);
for i=1:length(event)
  % test whether they match with the selected arrays
  if testvalue && isnumeric(value),         sel(i) = sel(i) && any(event(i).value == value);               end
  if testvalue && ischar(value),             sel(i) = sel(i) && any(strcmp(event(i).value,value));          end
  if testsample,        sel(i) = sel(i) && any(event(i).sample == sample);        end
  if testtimestamp,     sel(i) = sel(i) && any(event(i).timestamp == timestamp);  end
  if testoffset,        sel(i) = sel(i) && any(event(i).offset == offset);        end
  if testduration,      sel(i) = sel(i) && any(event(i).duration == duration);    end
  % test whether they lie within the specified range
  if testminsample,     sel(i) = sel(i) && (event(i).sample >= minsample);        end
  if testmaxsample,     sel(i) = sel(i) && (event(i).sample <= maxsample);        end
  if testminduration,   sel(i) = sel(i) && (event(i).duration >= minduration);    end
  if testmaxduration,   sel(i) = sel(i) && (event(i).duration <= maxduration);    end
  if testmintimestamp,  sel(i) = sel(i) && (event(i).timestamp >= mintimestamp);  end
  if testmaxtimestamp,  sel(i) = sel(i) && (event(i).timestamp <= maxtimestamp);  end
  if testminnumber,     sel(i) = sel(i) && (event(i).number >= minnumber);        end
  if testmaxnumber,     sel(i) = sel(i) && (event(i).number <= maxnumber);        end
  % this is potentially the slowest test, hence do it the last
  if testtype,          sel(i) = sel(i) && any(strcmp(event(i).type, type));      end
end

event = event(sel);
