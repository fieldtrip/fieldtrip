function event3 = align_presentation(event1, options1, event2, options2, hdr, feedback)

% ALIGN_PRESENTATION is a helper function to align events from a NBS Presentation log
% files to MEG/EEG triggers, or to a sequence of BOLD volumes.
%
% Use as
%    events3 = align_events(event1, options1, event2, options2)
% where
%   event1 = events from NBS Presentation log file
%   event2 = events from the MEG/EEG trigger channel
% or
%   event1 = events from NBS Presentation log file
%   event2 = events corresponding to each volume of the BOLD sequence
%
% The input "options1" and "options2" variables specify how the events should be
% mapped to each other. The output "events3" variable corresponds to the events from
% NBS Presentation log, but with the time aligned to the MEG/EEG dataset or to the
% BOLD volumes.
%
% See also DATA2BIDS, FT_READ_EVENT, FT_DEFINETRIAL

% determine whether the data corresponds to MEG/EEG or BOLD volumes
ismri = all(strcmp({event2.type}, 'volume'));

% make a copy of the presentation events
event3 = event1;

% select the corresponding events from the presentation file
sel1   = select_event(event1, options1.eventtype, options1.eventvalue);
event1 = event1(sel1);

% select the corresponding events from the trigger channel or BOLD volumes
sel2   = select_event(event2, options2.eventtype, options2.eventvalue);
event2 = event2(sel2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with potential mismatch between number of events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismri
  
  ft_info('%d volumes, %d presentation events', length(event2), length(event1));
  if length(event2)>length(event1)
    % this happens when the MRI scanner keeps running while presentation has already been stopped
    n = length(event2)-length(event1);
    ft_warning('discarding last %d volumes for realignment of events', n);
    event2 = event2(1:end-n);
  elseif length(event1)>length(event2)
    % this happens when DICOM volumes that are represented in the presentation log have been deleted from disk
    n = length(event1)-length(event2);
    switch options2.skip
      case 'first'
        ft_warning('discarding first %d presentation events for realignment of events', n);
        event1 = event1((n+1):end);
      case 'last'
        ft_warning('discarding last %d presentation events for realignment of events', n);
        event1 = event1(1:end-n);
      case 'none'
        ft_error('not enough volumes to match the presentation events');
    end % case
  end
  
else
  
  ft_info('%d triggers, %d presentation events', length(event2), length(event1));
  if length(event2)>length(event1)
    % don't know how to solve the situation of more triggers than presentation events
    ft_error('inconsistent number: %d triggers, %d presentation events', length(event2), length(event1));
  elseif length(event1)>length(event2)
    n = length(event1)-length(event2);
    % This could happen when due to acquisition problems there is more than one
    % dataset. If this is a known case, options2.skip can be used. Note that this
    % only works if there are exactly two datasets, not more.
    switch options1.skip
      case 'first'
        ft_warning('discarding first %d presentation events for realignment of events', n);
        event1 = event1((n+1):end);
      case 'last'
        ft_warning('discarding last %d presentation events for realignment of events', n);
        event1 = event1(1:end-n);
      case 'none'
        ft_error('not enough triggers to match the presentation events');
    end % case
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% align the events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% predict the sample number from the timestamp
model     = polyfit([event1.timestamp], [event2.sample], 1);
estimated = polyval(model, [event1.timestamp]);

if istrue(feedback)
  subplot(2,1,1)
  hold on
  % presentation timestamps are expressed in units of 0.1 miliseconds
  plot([event1.timestamp]/1e4, [event2.sample], 'b.')
  plot([event1.timestamp]/1e4, estimated, 'ro')
  xlabel('presentation time (s)')
  if ismri
    ylabel('MRI volumes')
  else
    ylabel('data samples')
  end
  legend({'observed', 'predicted'})
  
  subplot(2,1,2)
  plot([event1.timestamp]/1e4, ([event2.sample]-estimated)/hdr.Fs, 'g.')
  xlabel('presentation time (s)')
  ylabel('difference (s)')
end

% estimate the time in seconds of all presentation events
estimated = polyval(model, [event3.timestamp]);
estimated = round(1000*estimated)/1000; % round to three decimals
for i=1:numel(estimated)
  event3(i).sample = estimated(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sel = select_event(event, eventtype, eventvalue)
% this code is copied from FT_TRIALFUN_GENERAL

% start by selecting all events
sel = true(1, length(event)); % this should be a row vector

% select all events of the specified type
if ~isempty(eventtype)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).type, eventtype);
  end
end

% select all events with the specified value
if ~isempty(eventvalue)
  for i=1:numel(event)
    sel(i) = sel(i) && ismatch(event(i).value, eventvalue);
  end
end
