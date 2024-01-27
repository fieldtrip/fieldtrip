classdef MultiscaleElectrophysiologyData_1p0 < MultiscaleElectrophysiologyData
    % Class MULTISCALEELECTROPHYSIOLOGYDATA_1P0 process MED 1.0 channel data
    %
    % Syntax:
    %   this = MultiscaleElectrophysiologyDile_1p0();
    %   this = __(wholename);
    %   this = __(filepath, filename);
    %   this = __(__, 'Level1Password', level_1_pw);
    %   this = __(__, 'Level2Password', level_2_pw);
    %   this = __(__, 'AccessLevel', access_level);
    %
    % Input(s):
    %   wholename       - [char] (optional) session fullpath plus channel
    %                     name of MED file
    %   filepath        - [str] (optional) fullpath of session recorded in
    %                     MEF file
    %   filename        - [str] (optional) name of MEF channel file,
    %                     including ext
    %   level_1_pw      - [str] (para) password of level 1 (default = '')
    %   level_2_pw      - [str] (para) password of level 2 (default = '')
    %   access_level    - [str] (para) data decode level to be used
    %                     (default = 1)
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

    % Copyright 2023 Richard J. Cui. Created: Sun 02/12/2023 10:20:18.872 PM
    % $Revision: 0.6 $  $Date: at 07/22/2023 10:48:01.801 PM $
    %
    % Rocky Creek Dr. NE
    % Rochester, MN 55906, USA
    %
    % Email: richard.cui@utoronto.ca

    % =====================================================================
    % properties
    % =====================================================================
    % MED file info
    % -------------
    properties (SetAccess = protected)
        ChannelMetadata (1, 1) struct % channel metadata structure
        % .records: event records
        % .metadata: channel metadata
        % .contigua: contiguity information
    end % properties

    properties (SetAccess = protected, Hidden = true)
        Level1Password (1, :) char % level 1 password
        Level2Password (1, :) char % [str] level 2 password
        AccessLevel (1, 1) double {mustBeInteger, ...
                                       mustBeMember(AccessLevel, [1, 2])} = 1 % [num] access level of data
    end % properties

    % =====================================================================
    % methods
    % =====================================================================
    % the constructor
    % ----------------
    methods

        function this = MultiscaleElectrophysiologyData_1p0(file1st, file2nd, options)
            % constructor of class MultiscaleElectrophysiologyData_1p0
            % --------------------------------------------------------
            % * parse input
            arguments
                file1st (1, :) char = '';
                file2nd (1, :) char = '';
            end % positional

            arguments
                options.Level1Password (1, :) char = '';
                options.Level2Password (1, :) char = '';
                options.AccessLevel (1, 1) double {mustBeInteger, mustBePositive} = 1;
            end % optional

            if isempty(file1st)

                filepath = '';
                filename = '';
            else

                if isempty(file2nd)
                    [fp, fn, ext] = fileparts(file1st);
                    filepath = fp;
                    filename = [fn, ext];
                else
                    filepath = file1st;
                    filename = file2nd;
                end % if

                pw_1 = options.Level1Password;
                pw_2 = options.Level2Password;
                access_level = options.AccessLevel;

            end % if

            % operations during construction
            % ------------------------------
            % * call superclass constructor(s)
            this@MultiscaleElectrophysiologyData();

            % * set and check MED version
            if isnan(this.MEDVersion)
                this.MEDVersion = 1.0;
            elseif this.MEDVersion ~= 1.0
                error('MultiscaleElectrophysiologyData_1p0:InvalidMEDVer', ...
                    'Invalid MED version %f. Can only handle 1.0.', this.MEDVersion);
            end % if

            % * set channel info
            if ~isempty(filepath) && ~isempty(filename)
                this.FilePath = filepath;
                this.FileName = filename;
                this.Level1Password = pw_1;
                this.Level2Password = pw_2;
                this.AccessLevel = access_level;

                % read channel metadata
                switch this.AccessLevel
                    case 1
                        password = this.Level1Password;
                    case 2
                        password = this.Level2Password;
                end % switch

                wholename = fullfile(this.FilePath, this.FileName);

                if isfolder(wholename) == false
                    error('MultiscaleElectrophysiologyData_1p0:invalidChannel', ...
                        'cannot find channel %s', wholename)
                end % if

                this.read_channel_metadata(wholename, password);

                % set sampling information
                this.ChanSamplingFreq = this.ChannelMetadata.metadata.sampling_frequency;
                this.getSampleTimeInterval();

                % analyze continuity
                this.analyzeContinuity();

            end % if

        end

    end % methods

    % other methods
    % -------------
    methods
        channel = read_channel_metadata(this, varargin) % read channel metadata
        seg_cont = analyzeContinuity(this, varargin) % analyze continuity of data sampling
        [x, t] = importSignal(this, varargin) % input MEF 3.0 time series channel
        data = read_med_ts_data_1p0(this, varargin) % lower level of reading MED 1.0 time series data from one channel
        pw = processPassword(this, varargin) % process MEF 3.0 password
    end % methods

end % classdef

% [EOF]
