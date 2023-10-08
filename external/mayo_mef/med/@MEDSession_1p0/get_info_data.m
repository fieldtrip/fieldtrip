function [sess_info, unit] = get_info_data(this)
    % MEDSESSION_1P0.GET_INFO_DATA get session info from data
    %
    % Syntax:
    %   [sessionifo, unit] = get_info_data(this)
    %
    % Input(s):
    %   sess_info       - [table] N x 16 tabel: 'ChannelName', 'SamplingFreq',
    %                     'Begin', 'Stop', 'Samples' 'IndexEntry',
    %                     'DiscountinuityEntry', 'SubjectEncryption',
    %                     'SessionEncryption', 'DataEncryption', 'Version',
    %                     'Institution', 'SubjectID', 'AcquistitionSystem',
    %                     'CompressionAlgorithm', 'Continuity', where N is the
    %                     number of channels.
    %   unit            - [str] unit of begin_stop: 'Index' (default), 'uUTC',
    %                     'Second', 'Minute', 'Hour', and 'Day'
    %
    % Output(s):
    %
    % Example:
    %
    % Note:
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Tue 02/21/2023 11:27:28.805 PM
    % $Revision: 0.2 $  $Date: Fri 10/06/2023 10:46:53.182 PM $
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
    end % positional

    % =========================================================================
    % main
    % =========================================================================
    var_names = {'ChannelName', 'SamplingFreq', 'Begin', 'Stop', 'Samples', ...
                     'IndexEntry', 'DiscountinuityEntry', 'SubjectEncryption', ...
                     'SessionEncryption', 'DataEncryption', 'Version', 'Institution', ...
                     'SubjectID', 'AcquisitionSystem', 'CompressionAlgorithm', 'Continuity'};
    var_types = {'string', 'double', 'double', 'double', 'double', 'double', ...
                     'double', 'logical', 'logical', 'logical', 'string', 'string', 'string', ...
                     'string', 'string', 'cell'};

    % get the meta data
    % -----------------
    metadata = this.MetaData;

    if isempty(metadata) == true % no metadata
        metadata = this.read_med_session_metadata_1p0();
        this.MetaData = metadata;
    end % if

    % get session info
    % ----------------
    chan_names = metadata.channel_name;
    num_chan = length(chan_names);

    if num_chan < 1 % no time series channel
        sess_info = table;
        unit = '';
    else % if
        unit = 'uutc';
        sz = [num_chan, numel(var_names)]; % size of the table of session info
        fp = this.SessionPath; % session path of channels
        sess_info = table('size', sz, 'VariableTypes', var_types, ...
            'VariableNames', var_names);

        for k = 1:num_chan
            chan_name_k = chan_names(k);
            fn_k = chan_name_k + ".ticd";
            ch_k = MultiscaleElectrophysiologyData_1p0(fp, fn_k, ...
                'Level1Password', this.Password.Level1Password, ...
                'Level2Password', this.Password.Level2Password, ...
                'AccessLevel', this.Password.AccessLevel);
            info_k = ch_k.ChannelMetadata.metadata;
            seg_cont_k = ch_k.analyzeContinuity; % continuity table of the channel

            % assign values
            sess_info.ChannelName(k) = chan_name_k;
            sess_info.SamplingFreq(k) = info_k.sampling_frequency;
            sess_info.Begin(k) = info_k.start_time;
            sess_info.Stop(k) = info_k.end_time;
            sess_info.Samples(k) ...
                = info_k.absolute_end_sample_number - info_k.absolute_start_sample_number + 1; % number of samples in a segment
            % number of data blocks in each channel
            sess_info.IndexEntry(k) = nan;
            sess_info.DiscountinuityEntry(k) = height(seg_cont_k);

            % TODO: NA in MED 1.0
            sess_info.SubjectEncryption(k) = false;
            sess_info.SessionEncryption(k) = false;
            sess_info.DataEncryption(k) = false;

            sess_info.Version(k) = this.MEDVersion;
            sess_info.Institution(k) = info_k.recording_institution;
            sess_info.SubjectID(k) = info_k.subject_ID;

            % TODO: NA in MEF 3.0
            sess_info.AcquisitionSystem(k) = info_k.equipment_description;
            sess_info.CompressionAlgorithm(k) = 'not available';

            % continuity table
            sess_info.Continuity{k} = seg_cont_k;
        end % for

    end % if

end % function get_info_data

% [EOF]
