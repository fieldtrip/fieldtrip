function [x, t] = importSignal(this, start_end, st_unit, filepath, filename, options)
    % MULTISCALEELECTROPHYSIOLOGYDATA_1P0.IMPORTSIGNAL imports MED 1.0 signal into MATLAB
    %
    % Syntax:
    %   [x, t] = importSignal(this)
    %   [x, t] = importSignal(__, start_end)
    %   [x, t] = importSignal(__, start_end, st_unit)
    %   [x, t] = importSignal(__, start_end, st_unit, filepath)
    %   [x, t] = importSignal(__, start_end, st_unit, filepath, filename)
    %   [x, t] = importSignal(__, 'Level1Password', level_1_pw)
    %   [x, t] = importSignal(__, 'Level2Password', level_2_pw)
    %   [x, t] = importSignal(__, 'AccessLevel', access_level)
    %
    % Input(s):
    %   this            - [obj] MultiscaleElectrophysiologyData_1p0 object
    %   start_end       - [1 x 2 array] (opt) [start time/index, end time/index] of
    %                     the signal to be extracted fromt the file, relative time (default:
    %                     the entire signal)
    %   st_unit         - [str] (opt) unit of start_end: 'index' (default), 'uutc',
    %                     'second', 'minute', 'hour', and 'day'
    %   filepath        - [str] (opt) directory of the session
    %   filename        = [str] (opt) filename of the channel
    %   level_1_pw      - [str] (para) password of level 1 (default = this.Level1Password)
    %   level_2_pw      - [str] (para) password of level 2 (default = this.Level2Password)
    %   access_level    - [str] (para) data decode level to be used
    %                     (default = this.AccessLevel)
    %
    % Output(s):
    %   x               - [num array] extracted signal
    %   t               - [num array] time indices of the signal in the file
    %
    % Example:
    %
    % Note:
    %   Import data from one channel of MED 1.0 file into MatLab.
    %
    % References:
    %
    % See also .

    % Copyright 2023 Richard J. Cui. Created: Thu 04/20/2023 12:24:33.818 AM
    % $Revision: 0.9 $  $Date: Tue 10/10/2023 10:57:29.531 PM $
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
        start_end (1, 2) double {mustBeNonnegative} = [0, inf]
        st_unit (1, 1) string {mustBeMember(st_unit, ["index", "uutc", "second", "minute", "hour", "day"])} = "index"
        filepath (1, :) string = string.empty(1, 0)
        filename (1, :) string = string.empty(1, 0)
    end % positional

    arguments
        options.Level1Password (1, 1) string = this.Level1Password
        options.Level2Password (1, 1) string = this.Level2Password
        options.AccessLevel (1, 1) double = this.AccessLevel
    end % optional

    % ======================================================================
    % main
    % ======================================================================
    % parameters
    % ----------
    l1_pw = options.Level1Password;
    l2_pw = options.Level2Password;
    al = options.AccessLevel;

    % password
    % --------
    if l1_pw == ""
        l1_pw = this.Level1Password;
    else
        this.Level1Password = l1_pw;
    end % if

    if l2_pw == ""
        l2_pw = this.Level2Password;
    else
        this.Level2Password = l2_pw;
    end % if

    if isnan(al) || isempty(al)
        al = this.AccessLevel;
    else
        this.AccessLevel = al;
    end % if

    pw = this.processPassword('Level1Password', l1_pw, ...
        'Level2Password', l2_pw, ...
        'AccessLevel', al);

    % get the channel metadata if both filepath and filename are provided
    % -------------------------------------------------------------------
    if isempty(filepath)
        filepath = this.FilePath;
    else
        this.FilePath = filepath;
    end % if

    if isempty(filename)
        filename = this.FileName;
    else
        this.FileName = filename;
    end % if

    wholename = fullfile(filepath, filename);

    if ~isempty(filepath) && ~isempty(filename)
        this.read_channel_metadata(wholename, pw);
    end % if

    ch_meta = this.ChannelMetadata.metadata;

    % start and end time points
    % -------------------------
    switch lower(st_unit)
        case "index"
            se_index = start_end;
        otherwise
            se_index = this.SampleTime2Index(start_end, st_unit = st_unit);
    end % switch-case

    start_ind = ch_meta.absolute_start_sample_number;
    end_ind = ch_meta.absolute_end_sample_number;

    if isempty(start_end)
        se_index = [start_ind, end_ind];
    end % if

    % check
    % -----
    if isnan(se_index(1)) == true || se_index(1) < start_ind
        se_index(1) = start_ind;
        fprintf('Warning: MultiscaleElectrophysiologyData_1p0:ImportSignal:discardSample:Reqested data samples are not valid or before the recording are discarded\n')
    end % if

    if isnan(se_index(2)) || se_index(2) > end_ind
        se_index(2) = end_ind;
        warning('MultiscaleElectrophysiologyData_1p0:ImportSignal:discardSample', ...
        'Reqested data samples are not valid or after the recording are discarded')
    end % if

    % verbose
    % -------
    num_samples = diff(start_end) + 1;

    if num_samples > 2 ^ 20
        verbo = true;
    else
        verbo = false;
    end % if

    % load the data
    % -------------
    if verbo
        [~, thisChannel] = fileparts(wholename);
        fprintf("-->Loading "+thisChannel + " ...")
        clear thisChannel
    end % if

    x = this.read_med_ts_data_1p0(wholename, 'Password', pw, ...
        'range_type', 'samples', ...
        'begin', se_index(1), ...
        'stop', se_index(2));
    x = double(x(:)).'; % change to row vector
    % find the indices corresponding to physically collected data
    if nargout == 2
        t = se_index(1):se_index(2);
    end % if

    if verbo, fprintf('Done!\n'), end % if

end % function importSignal

% ==========================================================================
% subroutines
% ==========================================================================

% [EOF]
