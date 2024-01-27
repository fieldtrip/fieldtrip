function [sess_info, unit] = get_info_data(this)
% MEFSESSION_2P1.GET_INFO_DATA get sesion information from data
%
% Syntax:
%   [sessionifo, unit] = get_info_data(this)
%
% Input(s):
%   this            - [obj] MEFSession_2p1 object
%
% Output(s):
%   unit            - [str] unit of begin_stop: 'Index' (default), 'uUTC',
%                     'Second', 'Minute', 'Hour', and 'Day'
%   sess_info       - [table] N x 14 tabel: 'ChannelName', 'SamplingFreq',
%                     'Begin', 'Stop', 'Samples' 'IndexEntry',
%                     'DiscountinuityEntry', 'SubjectEncryption',
%                     'SessionEncryption', 'DataEncryption', 'Version',
%                     'Institution', 'SubjectID', 'AcquistitionSystem',
%                     'CompressionAlgorithm', 'Continuity', where N is the
%                     number of channels.
%
% Example:
%
% Note:
%   This function obtains information about the session from the data
%   directly.  Other information, such as session directory and password,
%   should be provided via MEFSession_2p1 object.
%
% References:
%
% See also MEFSession_2p1.

% Copyright 2020 Richard J. Cui. Created: Fri 01/03/2020  4:19:10.683 PM
% $ Revision: 0.2 $  $ Date: Fri 02/07/2020 11:34:16.078 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% main
% =========================================================================
sess_path = this.SessionPath;
pw = this.Password;
var_names = {'ChannelName', 'SamplingFreq', 'Begin', 'Stop', 'Samples',...
    'IndexEntry', 'DiscountinuityEntry', 'SubjectEncryption',...
    'SessionEncryption', 'DataEncryption', 'Version', 'Institution',...
    'SubjectID', 'AcquisitionSystem', 'CompressionAlgorithm', 'Continuity'};
var_types = {'string', 'double', 'double', 'double', 'double', 'double',...
    'double', 'logical', 'logical', 'logical', 'string', 'string', 'string',...
    'string', 'string', 'cell'};

chan_list = dir(fullfile(sess_path, '*.mef')); % assume all channel data in one dir
if isempty(chan_list)
    sess_info = table;
    unit = '';
else % if
    unit = 'uUTC';
    num_chan = numel(chan_list); % number of channels
    sz = [num_chan, numel(var_names)];
    sess_info = table('size', sz, 'VariableTypes', var_types,...
        'VariableNames', var_names);
    for k = 1:num_chan
        fp_k = chan_list(k).folder;
        fn_k = chan_list(k).name;
        header_k = this.readHeader(fullfile(fp_k, fn_k), pw.Subject);
        mef_ver = sprintf('%d.%d', header_k.header_version_major,...
            header_k.header_version_minor);
        ch2_k = MultiscaleElectrophysiologyFile_2p1(fp_k, fn_k,...
            'SubjectPassword', pw.Subject,...
            'SessionPassword', pw.Session,...
            'DataPassword', pw.Data);
        seg_cont_k = ch2_k.analyzeContinuity;
        
        sess_info.ChannelName(k)  = header_k.channel_name;
        sess_info.SamplingFreq(k) = header_k.sampling_frequency;
        sess_info.Begin(k)        = header_k.recording_start_time;
        sess_info.Stop(k)         = header_k.recording_end_time;
        sess_info.Samples(k)      = header_k.number_of_samples;
        sess_info.IndexEntry(k)   = header_k.number_of_index_entries;
        sess_info.DiscountinuityEntry(k) = header_k.number_of_discontinuity_entries;
        sess_info.SubjectEncryption(k)   = header_k.subject_encryption_used;
        sess_info.SessionEncryption(k)   = header_k.session_encryption_used;
        sess_info.DataEncryption(k)      = header_k.data_encryption_used;
        sess_info.Version(k)      = mef_ver;
        sess_info.Institution(k)  = header_k.institution;
        sess_info.SubjectID(k)    = header_k.subject_id;
        sess_info.AcquisitionSystem(k) = header_k.acquisition_system;
        sess_info.CompressionAlgorithm(k) = header_k.compression_algorithm;
        sess_info.Continuity{k} = seg_cont_k;
    end % for
end % if

end

% [EOF]