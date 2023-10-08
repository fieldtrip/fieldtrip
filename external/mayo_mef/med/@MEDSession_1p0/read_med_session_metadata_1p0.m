function metadata = read_med_session_metadata_1p0(this, sess_path, password)
    % MEDSESSION_1P0.READ_MED_SESSION_METADATA_1P0 (summary)
    %
    % Syntax:
    %   metadata = read_mef_session_metadata_3p0(this)
    %   metadata = __(__, sess_path)
    %   metadata = __(__, sess_path, password)
    %
    % Input(s):
    %   this            - [obj] MEFSession_3p0 object
    %   sess_path     	- [str] (opt) path (absolute or relative) to the MEF3
    %                     session folder (default = path in this)
    %   password        - [struct] password structure to the MEF3 data (default
    %                     = password in this)
    %
    % Output(s):
    %   metadata    	- [struct] structure containing session metadata,
    %                     channels metadata, segments metadata and records
    %
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Tue 02/28/2023 12:20:45.413 AM
    % $Revision: 0.1 $  $Date: Thu 08/31/2023 11:46:14.616 PM $
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
        sess_path (1, 1) string = this.SessionPath
        password (1, 1) struct = this.Password
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    % get session info
    % ----------------
    pw = this.processPassword(password);

    return_channels = true;
    return_contigua = true;
    return_records = true;

    sess_info = MED_session_stats(sess_path, return_channels, return_contigua, ...
        return_records, pw);

    % get metadata
    % ------------
    % get channel names
    channels = sess_info.channels;
    num_channels = length(channels);
    chan_names = strings(1, num_channels);

    for k = 1:num_channels
        chan_names(k) = channels(k).metadata.channel_name;
    end % for k

    % metadata
    metadata = sess_info.metadata;
    metadata.channel_name = chan_names;

end % function read_med_session_metadata_1p0

% [EOF]
