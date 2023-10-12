function [path_to_sess, sess_name, sess_ext] = get_sess_parts(this, sess_path)
    % MEDSESSION.GET_SESS_PARTS get the parts of the session path
    %
    % Syntax:
    %
    %   [path_to_sess, sess_name, sess_ext] = get_sess_parts(this)
    %   __ = __(__, sess_path)
    %
    % Input(s):
    %   this            - [obj] MEFSession object
    %   sess_path       - [str] (opt) Session path, including session name and
    %                     session extension (default: this.SessionPath)
    %
    % Output(s):
    %   path_to_sess    - [str] path up to the session, not include session
    %                     name and extension
    %   sess_name       - [str] session name, no extension
    %   sess_ext        - [str] session extension, including '.'
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Tue 02/21/2023 10:47:52.431 PM
    % $Revision: 0.1 $  $Date: Tue 02/21/2023 10:47:52.436 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) MEDSession
        sess_path (1, :) char = this.SessionPath
    end % positional

    % ======================================================================
    % main
    % ======================================================================
    [path_to_sess, sess_name, sess_ext] = fileparts(sess_path);

    % update session info
    % -------------------
    this.PathToSession = path_to_sess;
    this.SessionName = sess_name;
    this.SessionExt = sess_ext;

end % function get_sess_parts

% [EOF]
