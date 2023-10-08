function evt = getEvent(this, channames)
    % medfieldtrip_1p0.GETEVENT get MED 1,0 events for FieldTrip
    %
    % Syntax:
    %   evt = getEvent(this)
    %   evt = getEvent(this, channames)
    %
    % Input(s):
    %   this            - [obj] MEDFieldTrip_1p0 object
    %   channames       - [string] channel names
    %
    % Output(s):
    %	evt             - [struct] FieldTrip event structure
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Tue 10/03/2023 12:24:28.720 AM
    % $Revision: 0.1 $  $Date: Tue 10/03/2023 12:24:28.786 AM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) MEDFieldTrip_1p0
        channames (1, :) string = []
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    if ~isempty(channames)
        sess_chan = this.ChannelName; % all channels in the session

        if all(ismember(channames, sess_chan)) == false
            error('MEFFieldTrip_3p0:getHeader:invalidChannel', ...
            'invalid channel names')
        end % if

    end % if

    % construct the event structure
    % -----------------------------
    cont = this.SessionContinuity; % table of session continuity
    num_evt = size(cont, 1); % number of events
    evt = struct([]); % initialize the event structure

    for k = 1:num_evt
        % type
        evt(k).type = 'continuity recording';
        % sample
        evt(k).sample = cont.start_index(k);
        % value
        evt(k).value = [];
        % offset
        evt(k).offset = cont.end_index(k);
        % duration (in samples)
        evt(k).duration = cont.end_index(k) - cont.start_index(k) + 1;
        % timestamp of sample start
        evt(k).timestamp = cont.start_time(k);
    end % for

end % function getEvent

% [EOF]
