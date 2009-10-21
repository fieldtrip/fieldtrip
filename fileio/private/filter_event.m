function event = filter_event(event, varargin)

% FILTER_EVENT does what its name implies
%
% Use as
%   event = filter_event(event, ...)
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
% See also READ_EVEN, WRITE_EVENT

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: filter_event.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.8  2009/01/14 08:47:51  roboos
% moved to fileio
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.6  2007/12/20 19:04:12  roboos
% added filtering for minnumber and maxnumber (in line with read_neuralyns_nev)
%
% Revision 1.5  2007/12/19 08:28:50  roboos
% fixed bug for those filters which accept a numeric array (inserted any(..) around teh comparison)
% added teh missing code for filtering on timestamp value
%
% Revision 1.4  2007/08/21 16:56:03  chrhes
% fixed some logic (typo) bugs in relation to setting the testing flags
%
% Revision 1.3  2007/08/16 13:33:20  chrhes
% fixed some typos in the if statements implementing the event selection
%
% Revision 1.2  2007/07/30 12:16:17  roboos
% updated documentation
% implemented min and max timestamp
%
% Revision 1.1  2007/06/13 14:47:35  roboos
% moved filter_event from fieldtrip to fileio module
%
% Revision 1.1  2007/06/06 12:41:22  roboos
% new implementation
%

% get the optional input arguments
type         = keyval('type', varargin);
value        = keyval('value', varargin);
sample       = keyval('sample', varargin);
timestamp    = keyval('timestamp', varargin);
offset       = keyval('offset', varargin);
duration     = keyval('duration', varargin);

% the numeric fields can also be filtered on a range
minsample    = keyval('minsample', varargin);
maxsample    = keyval('maxsample', varargin);
minduration  = keyval('minduration', varargin);
maxduration  = keyval('maxduration', varargin);
mintimestamp = keyval('mintimestamp', varargin);
maxtimestamp = keyval('maxtimestamp', varargin);
minnumber    = keyval('minnumber', varargin);
maxnumber    = keyval('maxnumber', varargin);

if ~isempty(type)
  % this can be specified as string or as cell-array, convert to cell-array
  if ~iscell(type)
    type = {type};
  end
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

% apply the filters
sel = true(length(event),1);
for i=1:length(event)
  % test whether they match with the selected arrays
  if testvalue,    sel(i) = sel(i) && any(event(i).value == value); end
  if testsample,   sel(i) = sel(i) && any(event(i).sample == sample); end
  if testtimestamp,sel(i) = sel(i) && any(event(i).timestamp == timestamp); end
  if testoffset,   sel(i) = sel(i) && any(event(i).offset == offset); end
  if testduration, sel(i) = sel(i) && any(event(i).duration == duration); end
  % test whether they lie within the specified range
  if testminsample,   sel(i) = sel(i) && (event(i).sample >= minsample); end
  if testmaxsample,   sel(i) = sel(i) && (event(i).sample <= maxsample); end
  if testminduration, sel(i) = sel(i) && (event(i).duration >= minduration); end
  if testmaxduration, sel(i) = sel(i) && (event(i).duration <= maxduration); end
  if testmintimestamp, sel(i) = sel(i) && (event(i).timestamp >= mintimestamp); end
  if testmaxtimestamp, sel(i) = sel(i) && (event(i).timestamp <= maxtimestamp); end
  if testminnumber,   sel(i) = sel(i) && (event(i).number >= minnumber); end
  if testmaxnumber,   sel(i) = sel(i) && (event(i).number <= maxnumber); end
  % this is potentially the slowest test, hence do it the last
  if testtype,     sel(i) = sel(i) && any(strcmp(event(i).type, type)); end
end

event = event(sel);