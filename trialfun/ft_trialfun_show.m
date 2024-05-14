function [trl, event] = ft_trialfun_show(cfg)

% FT_TRIALFUN_SHOW will show a summary of the event information on screen. It will
% not return an actual trial definition. This function should in general not be
% called directly, it will be called by FT_DEFINETRIAL.
%
% Use this function by calling
%   [cfg] = ft_definetrial(cfg)
% where the configuration structure should contain
%   cfg.dataset   = string with the filename
%   cfg.trialfun  = 'ft_trialfun_show'
%
% See also FT_DEFINETRIAL, FT_TRIALFUN_GENERAL, FT_TRIALFUN_GUI

% Copyright (C) 2005-2021, Robert Oostenveld
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

% most defaults are in trialdef
cfg.trialdef = ft_getopt(cfg, 'trialdef', struct());

% specify the default file formats
cfg.eventformat   = ft_getopt(cfg, 'eventformat');
cfg.headerformat  = ft_getopt(cfg, 'headerformat');
cfg.dataformat    = ft_getopt(cfg, 'dataformat');

% construct the low-level options as key-value pairs, these are passed to FT_READ_EVENT
eventopt = {};
eventopt = ft_setopt(eventopt, 'headerformat',  ft_getopt(cfg, 'headerformat'));        % is passed to low-level function, empty implies autodetection
eventopt = ft_setopt(eventopt, 'dataformat',    ft_getopt(cfg, 'dataformat'));          % is passed to low-level function, empty implies autodetection
eventopt = ft_setopt(eventopt, 'eventformat',   ft_getopt(cfg, 'eventformat'));         % is passed to low-level function, empty implies autodetection
eventopt = ft_setopt(eventopt, 'readbids',      ft_getopt(cfg, 'readbids'));
eventopt = ft_setopt(eventopt, 'detectflank',   ft_getopt(cfg.trialdef, 'detectflank'));
eventopt = ft_setopt(eventopt, 'trigshift',     ft_getopt(cfg.trialdef, 'trigshift'));
eventopt = ft_setopt(eventopt, 'chanindx',      ft_getopt(cfg.trialdef, 'chanindx'));
eventopt = ft_setopt(eventopt, 'threshold',     ft_getopt(cfg.trialdef, 'threshold'));
eventopt = ft_setopt(eventopt, 'tolerance',     ft_getopt(cfg.trialdef, 'tolerance'));
eventopt = ft_setopt(eventopt, 'combinebinary', ft_getopt(cfg.trialdef, 'combinebinary'));

% define a message identifier and overrule the uses-specified default for ft_info
me = 'FieldTrip:ft_trialfun_show:showevents';
ft_info('on', me);

% get the events
if isfield(cfg, 'event')
  ft_info(me, 'using the events from the configuration structure\n');
  event = cfg.event;
else
  ft_info(me, 'reading the events from ''%s''\n', cfg.headerfile);
  event = ft_read_event(cfg.headerfile, eventopt{:});
end

if isempty(event)
  ft_info(me, 'no events were found in the data\n');
else
  ft_info(me, 'the following events were found in the data:\n\n');
  eventtype = unique({event.type});
  for i=1:length(eventtype)
    sel = find(strcmp(eventtype{i}, {event.type}));
    try
      value         = {event(sel(:)).value};
      eventvalue    = unique(value);                         % cell-array with string value
      eventvaluestr = sprintf('''%s'' ', eventvalue{:});     % translate into a single string
    catch
      value         = cell2mat({event(sel(:)).value});
      eventvalue    = unique(value);                         % array with numeric values or empty
      eventvaluestr = num2str(eventvalue);                   % translate into a single string
    end
    
    % count the number of occurrences per unique value
    n = zeros(numel(eventvalue), 1);
    nstr = '';
    for k = 1:numel(eventvalue)
      if iscell(eventvalue)
        n(k) = sum(strcmp(value, eventvalue{k}));
      else
        n(k) = sum(value==eventvalue(k));
      end
    end
    if isempty(eventvalue)
      n = numel(sel);
    end
    nstr = num2str(n(:)');  

    if ~iscell(eventvalue) && ~isempty(eventvalue)
      % try to make a nice formatting for the display, for cellarrays it's too difficult for now
      tmp = num2str([eventvalue(:)';n(:)'], '%d ');
      eventvaluestr = tmp(1,:);
      nstr = tmp(2,:);
    end
    
    ft_info(me, 'event type: ''%s'' ', eventtype{i});
    if ~isempty(eventvalue)
      ft_info(me, 'with event values : %s', eventvaluestr);
      ft_info(me, 'these values occur: %s times', nstr);
    else 
      ft_info(me, 'with event_values : %s', '[]');
      ft_info(me, 'these values occur: %s times', nstr);
    end
    ft_info(me, '\n');
  end
end

% this function always returns an empty trial definition
trl = [];
