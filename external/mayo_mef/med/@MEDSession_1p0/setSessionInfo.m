function metadata = setSessionInfo(this, sesspath, password, sort_channel)
    % MEDSESSION_1P0.SETSESSIONINFO Set session information of MEDSession_1p0
    %
    % Syntax:
    %   setSessionInfo(this, sesspath, password)
    %   setSessionInfo(__, sort_channel)
    %
    % Input(s):
    %   this            - [obj] MEDSession_1p0 object
    %   sesspath        - [char] session path of MED 3.0 data
    %   password        - [struct] MED 1.0 password structure (see MEDSession_1p0
    %                     for the details)
    %   sort_channel    - [char] (opt) sort channel according to either 'alphabet' of
    %                     the channel names or 'number' of the acquisiton
    %                     channel number (default = 'alphabet')
    %
    % Output(s):
    %   metadata        - [struct] MED 1.0 session metadata
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Tue 02/14/2023 11:36:09.573 PM
    % $Revision: 0.1 $  $Date: Tue 09/12/2023 11:03:56.620 PM $
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
        sesspath (1, :) char
        password (1, 1) struct
        sort_channel (1, :) logical = false
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    this.SessionPath = sesspath; % set session path of MED data
    this.Password = password; % set password of MED data
    this.get_sess_parts(); % get session parts
    this.get_sessinfo; % get session information

    % sorting channels
    % ----------------
    metadata = this.MetaData;
    ch_names = metadata.channel_name;

    if sort_channel == true
        ch_names = sort(ch_names);
        metadata.ChannelName = ch_names;
    end % if

    this.ChannelName = ch_names;
    this.MetaData = metadata;
end % function setSessionInfo

% [EOF]
