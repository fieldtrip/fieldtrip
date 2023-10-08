function varargout = get_sessinfo(this)
    % MEDSESSION.GET_SESSINFO get session information from MED data
    %
    % Syntax:
    %   [channame, start_end, unit, sess_info] = get_sessinfo(this)
    %
    % Input(s):
    %   this            - [obj] MEDSession object
    %
    % Output(s):
    %   channame        - [string] the name(s) of the data channel in the
    %                     directory of session.
    %   begin_stop      - [1 x 2 array] [begin time/index, stop time/index] of
    %                     the entire signal
    %   unit            - [str] unit of begin_stop: 'Index' (default), 'uUTC',
    %                     'Second', 'Minute', 'Hour', and 'Day'
    %   sess_info       - [table] N x 13 tabel: 'ChannelName', 'SamplingFreq',
    %                     'Begin', 'Stop', 'Samples' 'IndexEntry',
    %                     'DiscountinuityEntry', 'SubjectEncryption',
    %                     'SessionEncryption', 'DataEncryption', 'Version',
    %                     'Institution', 'SubjectID', 'AcquistitionSystem',
    %                     'CompressionAlgorithm', where N is the number of
    %                     channels.
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also MEDSession_1p0, get_info_data.

    % Copyright 2023 Richard J. Cui. Created: Tue 02/21/2023 11:21:04.215 PM
    % $Revision: 0.1 $  $Date: Mon 09/11/2023 10:30:51.387 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % Main process
    % ======================================================================
    % table of channel info
    % ---------------------
    [sess_info, unit] = this.get_info_data;

    if isempty(sess_info)
        si_value = get_sess_info_value(sess_info, unit, true);
        cprintf([1, 0.5, 0], 'The session is likely empty.\n')
    else
        this.SessionInformation = sess_info;
        si_value = get_sess_info_value(sess_info, unit, false);
    end % if

    % update paras of MEDSession
    % --------------------------
    this.ChannelName = si_value.channame;
    this.SamplingFrequency = si_value.fs;
    this.Samples = si_value.samples;
    this.DataBlocks = si_value.num_data_block;
    this.TimeGaps = si_value.num_time_gap;
    this.BeginStop = si_value.begin_stop;
    this.Unit = si_value.unit;
    this.Institution = si_value.institution;
    this.SubjectID = si_value.subj_id;
    this.AcquisitionSystem = si_value.acq_sys;
    this.CompressionAlgorithm = si_value.comp_alg;
    this.SessionInformation = si_value.sess_info;
    this.SessionContinuity = si_value.cont;

    % ======================================================================
    % Output
    % ======================================================================
    if nargout > 0
        varargout{1} = channame;
    end % if

    if nargout > 1
        varargout{2} = begin_stop;
    end % if

    if nargout > 2
        varargout{3} = unit;
    end % if

    if nargout > 3
        varargout{4} = sess_info;
    end % if

end % function get_sessinfo

% ==========================================================================
% Subroutines
% ==========================================================================
function si_value = get_sess_info_value(sess_info, unit, set_empty)

    if set_empty == true
        si_value = struct('channame', '', ...
            'fs', NaN, ...
            'samples', NaN, ...
            'num_data_block', [], ...
            'begin_stop', [], ...
            'unit', '', ...
            'institution', '', ...
            'subj_id', '', ...
            'acq_sys', '', ...
            'comp_alg', '', ...
            'sess_info', struct([]), ...
            'cont', table([]));
    else
        si_value.channame = sess_info.ChannelName';
        si_value.fs = unique(sess_info.SamplingFreq);
        si_value.samples = unique(sess_info.Samples);
        si_value.num_data_block = unique(sess_info.IndexEntry);
        si_value.num_time_gap = unique(sess_info.DiscountinuityEntry) - 1;
        si_value.begin_stop = [unique(sess_info.Begin), unique(sess_info.Stop)];
        si_value.unit = unit;
        si_value.institution = unique(sess_info.Institution);
        si_value.subj_id = unique(sess_info.SubjectID);
        si_value.acq_sys = unique(sess_info.AcquisitionSystem);
        si_value.comp_alg = unique(sess_info.CompressionAlgorithm);
        si_value.sess_info = sess_info;
        si_value.cont = sess_info.Continuity{1}; % use the one from 1st channel
    end % if

end % function

% [EOF]
