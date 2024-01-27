function [X, t] = import_sess(this, begin_stop, bs_unit, sel_chan, options)
    % MEDSESSION_1P0.IMPORT_SESS import session of MED 1.0 data
    %
    % Syntax:
    %   [X, t] = import_sess(this, begin_stop, bs_unit, sel_chan)
    %   [X, t] = import_sess(__, pw)
    %
    % Input(s):
    %   this            - [obj] MEDSession_1p0 object
    %   begin_stop      - [num] 1 x 2 array of begin and stop points of
    %                     importing the session, relative time points
    %   bs_unit         - [str] unit of begin_stop: 'uUTC','Index', 'Second',
    %                     'Minute', 'Hour', and 'Day'.
    %   sel_chan        - [str array] the names of the selected channels
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
    % Example:
    %
    % Note:
    %   Import data from different channels of the session.
    %
    % References:
    %
    % See also importSignal, importSession.

    % Copyright 2023 Richard J. Cui. Created: Thu 09/14/2023 11:56:16.600 PM
    % $Revision: 0.3 $  $Date: Tue 10/10/2023 10:57:29.531 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) MEDSession_1p0
        begin_stop (1, 2) double
        bs_unit (1, :) char ...
            {mustBeMember(bs_unit, {'uutc', 'index', 'second', 'minute', 'hour', 'day'})}
        sel_chan (1, :) string
    end % positional

    arguments
        options.Password (1, 1) struct = struct()
    end % name-value

    pw = options.Password;

    if isempty(pw)
        pw = this.Password;
    end % if

    sess_path = this.SessionPath;

    % ======================================================================
    % main
    % ======================================================================
    num_chan = numel(sel_chan);
    X = [];

    for k = 1:num_chan
        fn_k = sel_chan(k) + ".ticd";
        [x_k, t] = this.importSignal(begin_stop, bs_unit, sess_path, fn_k, ...
            'Level1Password', pw.Level1Password, ...
            'Level2Password', pw.Level2Password, ...
            'AccessLevel', pw.AccessLevel);
        x_k = x_k(:).'; % make sure it is a horizontal vector

        X = cat(1, X, x_k);
    end % for

end % function import_sess

% [EOF]
