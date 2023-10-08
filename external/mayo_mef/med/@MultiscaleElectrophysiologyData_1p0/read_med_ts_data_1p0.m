function data = read_med_ts_data_1p0(this, channel_path, options)
    % MULTISCALEELECTROPHYSIOLOGYDATA_1P0.READ_MED_TS_DATA_1P0 read the MED 1.0 data from a time-series channel
    %
    % Syntax:
    %
    % Input(s):
    %
    % Output(s):
    %
    % Example:
    %
    % Note:
    %   If the range_type is "time", then the begin and stop are in offset uUTC.
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Sat 07/22/2023 10:51:35.758 PM
    % $Revision: 0.2 $  $Date: Mon 08/21/2023 12:33:26.391 AM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % ======================================================================
    % parse inputs
    % ======================================================================
    arguments
        this (1, 1) MultiscaleElectrophysiologyData_1p0
        channel_path (1, 1) string {mustBeNonzeroLengthText} % the path to the channel
    end % positional

    arguments
        options.Password (1, :) char = '' % the password to the file
        options.range_type (1, 1) string {mustBeMember(options.range_type, ["time", "samples"])} = "samples" % the range type
        options.begin (1, 1) double
        options.stop (1, 1) double
    end % optional

    pw = options.Password;
    range_type = options.range_type;
    begin = options.begin;
    stop = options.stop;

    % ======================================================================
    % main
    % ======================================================================
    start_time = this.ChannelMetadata.metadata.start_time;

    if range_type == "time"
        % convert the begin and stop to index
        begin_index = this.SampleTime2Index(begin - start_time, st_unit = "uUTC");
        stop_index = this.SampleTime2Index(stop - start_time, st_unit = "uUTC");
    else
        begin_index = begin;
        stop_index = stop;
    end % if

    sess = read_MED(channel_path, [], [], begin_index, stop_index, pw);
    data = sess.channels.data;

end % function read_med_ts_data_1p0

% [EOF]
