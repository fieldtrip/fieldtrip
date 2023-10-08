function [X, t] = importSession(this, start_end, se_unit, sess_path, options)
    % MEDSESSION.importSession import MED session data
    %
    % Syntax:
    %   [X, t] = importSession(this)
    %   [X, t] = importSession(__, start_end)
    %   [X, t] = importSession(__, start_end, se_unit)
    %   [X, t] = importSession(__, start_end, se_unit, sess_path)
    %   [X, t] = importSession(__, 'SelectedChannel', sel_chan, 'Password', pw)
    %
    % Imput(s):
    %   this            - [obj] MEDSession object
    %   start_end       - [num] (opt) 1 x 2 array of begin and stop points of
    %                     importing the session, relative time points (default:
    %                     if it is not given, the system will try to use the
    %                     property this.StartEnd. if StartEnd is empty, the system
    %                     will read the entire signal)
    %   se_unit         - [str] (opt) unit of start_end: 'uUTC' (default),
    %                     'Index', 'Second', 'Minute', 'Hour', and 'Day'.
    %   sess_path       - [str] (opt) session path (default: this.SessionPath)
    %   sel_chan        - [str array] (para) the names of the selected channels
    %                     (default: this.SelectedChannel, otherwise all channels)
    %   pw              - [struct] (para) password structure
    %                     .Session      : session password
    %                     .Subject      : subject password
    %                     .Data         : data password
    %
    % Output(s):
    %   X               - [num array] M x N array, where M is the number of
    %                     channels and N is the number of signals extracted
    %   t               - [num] 1 x N array, time indeces of the signals
    %
    % Note:
    %   Import data from different channels of the session.
    %
    % See also importSignal.

    % Copyright 2020-2023 Richard J. Cui. Created: Wed 01/08/2020 11:16:21.943 PM
    % $Revision: 0.3 $  $Date: Fri 10/06/2023 12:15:37.649 AM $
    %
    % 1026 Rocky Creek Dr NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % =========================================================================
    % parse inputs
    % =========================================================================
    arguments
        this (1, 1) MEDSession
        start_end double = []
        se_unit char {mustBeMember(se_unit, ...
                          {'uutc', 'index', 'second', 'minute', 'hour', 'day'})} = 'uutc'
        sess_path char = []
    end % positional

    arguments
        options.SelectedChannel string = []
        options.Password struct = struct([])
    end % optional

    sel_chan = options.SelectedChannel;
    pw = options.Password;

    if isempty(start_end)

        if isempty(this.StartEnd)
            start_end = this.abs2relativeTimePoint(this.BeginStop, this.Unit); % absolute time points
            this.StartEnd = start_end;
            this.SEUnit = this.Unit;
        else
            start_end = this.StartEnd;
            se_unit = this.SEUnit;
        end % if

    else
        this.StartEnd = start_end;
        this.SEUnit = se_unit;
    end % if

    if isempty(sess_path)
        sess_path = this.SessionPath;

        if isempty(sess_path)
            warning('MEDSession:importSession:noSession', ...
            'No session is selected')
            X = [];
            t = [];
            return
        end % if

    else
        this.SessionPath = sess_path;
    end

    if isempty(sel_chan)
        sel_chan = this.SelectedChannel;

        if isempty(sel_chan)
            sel_chan = this.ChannelName;

            if isempty(sel_chan)
                warning('MEDSession:importSession:emptySession', ...
                'Either the session is empty or no channel has been selected')
                X = [];
                t = [];
                return
            end % if

            this.SelectedChannel = sel_chan;
        end % if

    else
        this.SelectedChannel = sel_chan;
    end % if

    if isempty(pw)
        pw = this.Password;
    end % if

    % ======================================================================
    % input session
    % ======================================================================
    begin_stop = this.relative2absTimePoint(start_end, se_unit); % to absolute time points
    [X, t] = this.import_sess(begin_stop, se_unit, sel_chan, Password = pw);

end % function

% =========================================================================
% subroutines
% =========================================================================

% [EOF]
