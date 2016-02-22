function [ev] = ft_recodeevent(cfg, event, trl)

% FT_RECODEEVENT will recode the event structure, given the trial
% definition that was analyzed
%
% In FieldTrip, you always start with defining a "trl" field containing
% the samples in the raw datafile that you want to analyze. That "trl"
% is based on the events in the dataset. After artifact rejection, it may
% be the case that trials have been removed completely, or that trials
% have been cut into pieces. This complicates finding a match between the
% original events and the pieces of data that are analyzed. This functino
% restores that match.
%
% Use as
%   [ev] = ft_recodeevent(cfg, data)
% where cfg is a structure with configuration settings and data contains the
% (nested) configuration that describes the original trial definition and
% event structure.
%
% Alternatively, you can also specify the event structure and trial definition
% yourself with
%   [ev] = ft_recodeevent(cfg, event, trl)
%
% the configuration can contain
%   cfg.eventtype   = empty, 'string' or cell-array with multiple strings
%   cfg.eventvalue  = empty or a list of event values (can be numeric or string)
%
%   cfg.searchrange = 'anywhere'      search anywhere for the event, (default)
%                     'insidetrial'   only search inside
%                     'outsidetrial'  only search outside
%                     'beforetrial'   only search before the trial
%                     'aftertrial'    only search after  the trial
%                     'beforezero'    only search before time t=0 of each trial
%                     'afterzero'     only search after  time t=0 of each trial
%
%   cfg.nearestto   = 'trialzero'     compare with time t=0 for each trial (default)
%                     'trialbegin'    compare with the begin of each trial
%                     'trialend'      compare with the end of each trial
%
%   cfg.match       = 'exact' or 'nearest'
%
%   cfg.output      = 'event'             the event itself
%                     'eventvalue'        the value of the event
%                     'eventnumber'       the number of the event
%                     'samplenumber'      the sample at which the event is located
%                     'samplefromoffset'  number of samples from t=0 (c.f. response time)
%                     'samplefrombegin'   number of samples from the begin of the trial
%                     'samplefromend'     number of samples from the end   of the trial
%
% See also FT_DEFINETRIAL, FT_REDEFINETRIAL, FT_PREPROCESSING

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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.eventtype   = ft_getopt(cfg, 'eventtype',   []);
cfg.eventvalue  = ft_getopt(cfg, 'eventvalue',  []);
cfg.searchrange = ft_getopt(cfg, 'searchrange', 'anywhere');
cfg.nearestto   = ft_getopt(cfg, 'nearestto',   'trialzero');
cfg.match       = ft_getopt(cfg, 'match',       'nearest');
cfg.output      = ft_getopt(cfg, 'output',      'eventvalue');

% these should be numeric lists or cell-arrays with strings
if ischar(cfg.eventtype)
  cfg.eventtype = {cfg.eventtype};
end
if ischar(cfg.eventvalue)
  cfg.eventvalue = {cfg.eventvalue};
end

if nargin==2
  % event and trl are not specified in the function call, but the data is given ->
  % try to locate event and trl in the configuration
  data  = event;                       % rename the input variable
  if isfield(data, 'cfg')
    event = ft_findcfg(data.cfg, 'event');  % search for the event field
    trl   = ft_findcfg(data.cfg, 'trl');    % search for the trl field
  else
    event = [];
    trl   = [];
  end
  if isempty(event)
    error('could not locate event structure in the data');
  elseif isempty(trl)
    error('could not locate trial definition in the data');
  end
elseif nargin~=3
  error('incorrect number of input arguments');
end

Ntrl   = size(trl,1);
Nevent = length(event);

% select the events of interest
fprintf('trial definition describes %d trials\n', Ntrl);
fprintf('original event structure contains %d events\n', Nevent);
selecttype  = zeros(Nevent,1);
selectvalue = zeros(Nevent,1);
for i=1:Nevent
  % test whether this event should be selected
  if ~isempty(cfg.eventtype)
    selecttype(i) = ~isempty(intersect(cfg.eventtype, event(i).type));
  else
    selecttype(i) = 1;
  end
  % test whether this event should be selected
  if ~isempty(cfg.eventvalue)
    selectvalue(i) = ~isempty(intersect(cfg.eventvalue, event(i).value));
  else
    selectvalue(i) = 1;
  end
end
fprintf('selected %d events based on event type\n', sum(selecttype));
fprintf('selected %d events based on event value\n', sum(selectvalue));
fprintf('selected %d events based on event type and value\n', sum(selecttype.*selectvalue));
eventnum = find(selecttype.*selectvalue);
event = event(eventnum);
Nevent = length(event);

if Nevent<1
  error('there are no events to analyze');
end

% make a list with the sample, offset and duration of each event
% and sort the events according to the sample at which they occurred
sample   = zeros(Nevent,1);
offset   = zeros(Nevent,1);
duration = zeros(Nevent,1);
for i=1:Nevent
  sample(i) = event(i).sample;
  if ~isempty(event(i).offset)
    offset(i) = event(i).offset;
  else
    offset(i) = nan;
  end
  if ~isempty(event(i).duration)
    duration(i) = event(i).duration;
  else
    duration(i) = nan;
  end
end
[sample, indx] = sort(sample);   % sort the samples
offset = offset(indx);           % sort the offset accordingly
duration = duration(indx);       % sort the duration accordingly
event = event(indx);             % sort the events accordingly
eventnum = eventnum(indx);       % sort the numbers of the original events

for i=1:Ntrl
  trlbeg    = trl(i,1);
  trlend    = trl(i,2);
  trloffset = trl(i,3);
  trlzero   = trlbeg - trloffset;  % the sample that corresponds with t=0

  if strcmp(cfg.nearestto, 'trialzero')
    trlsample = trlzero;             % the sample that corresponds with t=0
  elseif strcmp(cfg.nearestto, 'trialbegin')
    trlsample = trlbeg;              % the sample at which the trial begins
  elseif strcmp(cfg.nearestto, 'trialend')
    trlsample = trlend;              % the sample at which the trial ends
  else
    error('incorrect specification of cfg.nearestto')
  end

  % compute a "distance" measure for each event towards this trial
  switch cfg.searchrange
    case 'anywhere'
      distance = abs(sample - trlsample);
    case 'beforezero'
      distance = abs(sample - trlsample);
      distance(find(sample>=trlzero)) = inf;
    case 'afterzero'
      distance = abs(sample - trlsample);
      distance(find(sample<=trlzero)) = inf;
    case 'beforetrial'
      distance = abs(sample - trlsample);
      distance(find(sample>=trlbeg)) = inf;
    case 'aftertrial'
      distance = abs(sample - trlsample);
      distance(find(sample<=trlend)) = inf;
    case 'insidetrial'
      distance = abs(sample - trlsample);
      distance(find((sample<trlbeg) | (sample>trlend))) = inf;
    case 'outsidetrial'
      distance = abs(sample - trlsample);
      distance(find((sample>=trlbeg) & (sample<=trlend))) = inf;
    otherwise
      error('incorrect specification of cfg.searchrange');
  end

  % determine the event that has the shortest distance towards this trial
  [mindist, minindx] = min(distance);
  if length(find(distance==mindist))>1
    error('multiple events are at the same distance from the trial');
  end

  if isinf(mindist)
    % no event was found
    ev(i) = nan;
  elseif mindist~=0 && strcmp(cfg.match, 'exact')
    % the event is not an exact match
    ev(i) = nan;
  else
    switch cfg.output
      case 'event'
        ev(i) = event(minindx);
      case 'eventvalue'
        ev(i) = event(minindx).value;
      case 'eventnumber'
        ev(i) = eventnum(minindx);
      case 'samplenumber'
        ev(i) = event(minindx).sample;
      case 'samplefromoffset'
        ev(i) = event(minindx).sample - trlzero;
      case 'samplefrombegin'
        ev(i) = event(minindx).sample - trlbeg;
      case 'samplefromend'
        ev(i) = event(minindx).sample - trlend;
      otherwise
        error('incorrect specification of cfg.output');
    end
  end

end % looping over all trials
