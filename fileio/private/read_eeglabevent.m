% read_eeglabevent() - import EEGLAB dataset events
%
% Usage:
%     >> event = read_eeglabevent(filename, ...);
%
% Inputs:
%     filename - [string] file name
%
% Optional inputs:
%     'header' - FILEIO structure header
%
% Outputs:
%     event - FILEIO toolbox event structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2008-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

function event = read_eeglabevent(filename, varargin)

if nargin < 1
  help read_eeglabheader;
  return;
end;

hdr = ft_getopt(varargin, 'header');

if isempty(hdr)
  hdr = read_eeglabheader(filename);
end

event    = [];                % these will be the output in FieldTrip format
oldevent = hdr.orig.event;    % these are in EEGLAB format

missingFieldFlag=false;

if ~isempty(oldevent)
    if ~isfield(oldevent,'code') && ~isfield(oldevent,'value')  && ~isfield(oldevent,'setname')
        disp('Warning: No ''value'' field in the events structure.');
        missingFieldFlag=true;
    end;
    
    if ~isfield(oldevent,'type')
        disp('Warning: No ''type'' field in the events structure.');
        missingFieldFlag=true;
    end;
    
    if missingFieldFlag
        if ~isfield(oldevent,'setname') %accommodate Widmann's pop_grandaverage function
            disp('EEGlab data files should have both a ''value'' field');
            disp('to denote the generic type of event, as in ''trigger'', and a ''type'' field');
            disp('to denote the nature of this generic event, as in the condition of the experiment.');
            disp('Note also that this is the reverse of the FieldTrip convention.');
        end;
    end;
end;

for index = 1:length(oldevent)

  if isfield(oldevent,'code')
    type = oldevent(index).code;
  elseif isfield(oldevent,'value')
    type = oldevent(index).value;
  else
    type = 'trigger';
  end;

  % events can have a numeric or a string value
  if isfield(oldevent,'type')
    value  = oldevent(index).type;
  else
    value = 'default';
  end;

  % this is the sample number of the concatenated data to which the event corresponds
  sample = oldevent(index).latency;

  % a non-zero offset only applies to trial-events, i.e. in case the data is
  % segmented and each data segment needs to be represented as event. In
  % that case the offset corresponds to the baseline duration (times -1).
  offset = 0;

  if isfield(oldevent, 'duration')
    duration = oldevent(index).duration;
  else
    duration = 0;
  end;

  % add the current event in fieldtrip format
  event(index).type     = type;     % this is usually a string, e.g. 'trigger' or 'trial'
  event(index).value    = value;    % in case of a trigger, this is the value
  event(index).sample   = sample;   % this is the sample in the datafile at which the event happens
  event(index).offset   = offset;   % some events should be represented with a shifted time-axix, e.g. a trial with a baseline period
  event(index).duration = duration; % some events have a duration, such as a trial

end;

if hdr.nTrials>1
  % add the trials to the event structure
  for i=1:hdr.nTrials
    event(end+1).type     = 'trial';
    event(end  ).sample   = (i-1)*hdr.nSamples + 1;
    if isfield(oldevent,'setname') && (length(oldevent) == hdr.nTrials)
        event(end  ).value    = oldevent(i).setname; %accommodate Widmann's pop_grandaverage function
    else
        event(end  ).value    = [];
    end;
    event(end  ).offset   = -hdr.nSamplesPre;
    event(end  ).duration =  hdr.nSamples;
  end
end

