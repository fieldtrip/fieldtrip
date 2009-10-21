function [event] = read_nmc_archive_k_event(eventfile)

% READ_NMC_ARCHIVE_K_EVENT extracts event-data from nmc_archive_k datasets
%
% Use as
%   event = read_nmc_archive_k_event(eventfile)
%
%
% This function specifically only reads data from one of the archived
% datasets of the Neurophysiological Mechanisms of Cognition group of
% Eric Maris, at the Donders Centre for Cognition, Radboud University,
% Nijmegen, the Netherlands. It should not be used for any other data
% format.
%
%
% First version: 2009/09/31 - roevdmei
%


% Checking events-file:
if exist(eventfile,'file') ~= 2
    error('no events.mat file found in specified directory');
end

% Load event-file as events_old
events_old = load(eventfile);

% Produce variables from eventfile-name
slashpos = strfind(eventfile, '/');
subjname = eventfile((slashpos(end-2)+1):(slashpos(end-1)-1));
sessionname = eventfile((slashpos(end)+1):end-10);

% Changing events to FieldTrip format and starting with empty event structure and with numbering events
event = [];
eventnum = 1;
for ievent = 1:length(events_old.events)
    slashpos = strfind(events_old.events(ievent).eegfile, '/');
    if ~isempty(events_old.events(ievent).eegfile) % dont select events that are not labelled to a specific session
        if  strcmp(events_old.events(ievent).eegfile((slashpos(end)+1):end),sessionname) % checking if session corresponds to requested session
            if ~isempty(events_old.events(ievent).eegoffset) % dont select events that contain no sample number
                if ~isfield(events_old.events(ievent), 'artifacts') || (isfield(events_old.events(ievent), 'artifacts') && isempty(events_old.events(ievent).artifacts)) % dont select events that contain predefined artifacts
                    if ~isfield(events_old.events(ievent), 'ispractice') || (isfield(events_old.events(ievent), 'ispractice') && ~isempty(events_old.events(ievent).ispractice) && events_old.events(ievent).ispractice == 0) % aviods using practice trials in analyses
                        event(eventnum).type        = events_old.events(ievent).type;
                        event(eventnum).sample      = events_old.events(ievent).eegoffset;
                        if iscell(events_old.events(ievent).mode)
                            event(eventnum).mode    = events_old.events(ievent).mode{1}; % .mode variables are a cell array in some subjects
                        else
                            event(eventnum).mode    = events_old.events(ievent).mode;
                        end
                        event(eventnum).item        = events_old.events(ievent).item;
                        event(eventnum).cueitem     = events_old.events(ievent).cueitem;
                        event(eventnum).correctresp = events_old.events(ievent).correct;
                        event(eventnum).istarget    = events_old.events(ievent).istarget;
                        event(eventnum).rt          = events_old.events(ievent).rt;
                        event(eventnum).probe_pos   = events_old.events(ievent).probe_pos;
                        event(eventnum).listlength  = events_old.events(ievent).listlen;
                        eventnum = eventnum + 1;
                    end % exist ispractice
                end % exist artifacts
            end % isempty
        end % strcmp
    end % isempty
end % ievent


% Send warning if no events are found for specified session
if isempty(event)
    warning(['no events found for session: ' sessionname ' of subject: ' subjname])
end