% mff_fileio_read_event() - import MFF events
%
% Usage:
%     >> event = mff_fileio_read_event(filename, ...);
%
% Inputs:
%     filename - [string] file name
%
% Optional inputs:
%     'header' - FILEIO structure header
%
% Outputs:
%     event - FILEIO toolbox event structure

% This file is part of mffmatlabio.
%
% mffmatlabio is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% mffmatlabio is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with mffmatlabio.  If not, see <https://www.gnu.org/licenses/>.

function event = mff_fileio_read_event(filename, varargin)

if nargin < 1
    help mff_fileio_read_event;
    return;
end

hdr = ft_getopt(varargin, 'header');

if isempty(hdr)
    hdr = mff_fileio_read_header(filename);
end

event    = [];                % these will be the output in FieldTrip format
oldevent = hdr.orig.event;    % these are in EEGLAB format

if ~isempty(oldevent)
    nameList=fieldnames(oldevent);
else
    nameList=[];
end
nameList=setdiff(nameList,{'type','value','sample','offset','duration','latency'});

for index = 1:length(oldevent)
    
    if isfield(oldevent,'code')
        type = oldevent(index).code;
    elseif isfield(oldevent,'value')
        type = oldevent(index).value;
    else
        type = 'trigger';
    end
    
    % events can have a numeric or a string value
    if isfield(oldevent,'type')
        value  = oldevent(index).type;
    else
        value = 'default';
    end
    
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
    end
    
    % add the current event in FieldTrip format
    event(index).type     = type;     % this is usually a string, e.g. 'trigger' or 'trial'
    event(index).value    = value;    % in case of a trigger, this is the value
    event(index).sample   = sample;   % this is the sample in the datafile at which the event happens
    event(index).offset   = offset;   % some events should be represented with a shifted time-axix, e.g. a trial with a baseline period
    event(index).duration = duration; % some events have a duration, such as a trial
    
    %add custom fields
    for iField=1:length(nameList)
        eval(['event(index).' nameList{iField} '=oldevent(index).' nameList{iField} ';']);
    end
    
end
% 
% if hdr.nTrials>1
%     % add the trials to the event structure
%     for i=1:hdr.nTrials
%         event(end+1).type     = 'trial';
%         event(end  ).sample   = (i-1)*hdr.nSamples + 1;
%         if isfield(oldevent,'setname') && (length(oldevent) == hdr.nTrials)
%             event(end  ).value    = oldevent(i).setname; %accommodate Widmann's pop_grandaverage function
%         else
%             event(end  ).value    = [];
%         end
%         event(end  ).offset   = -hdr.nSamplesPre;
%         event(end  ).duration =  hdr.nSamples;
%     end
% end
% 
