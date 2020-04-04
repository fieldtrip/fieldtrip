classdef MEFSession_3p0 < MEFSession & MultiscaleElectrophysiologyFile_3p0
    % Class MEFSESSION_3P0 process MEF 3.0 session
    % 
    % Syntax:
    %   this = MEFSession_3p0
    %   this = __(filename)
    %   this = __(filename, password)
    %
    % Input(s):
    %   filename    - [str] (opt) MEF 3.0 session path, channel or data file
    %                 (default = '')
    %   password    - [struct] (opt) structure of MEF 3.0 passowrd (default
    %                 = struct('Level1Password', '', 'Level2Password', '',...
    %                 'AccessLevel', 1);)
    %                 .Level1Password (default = '')
    %                 .Level2Password (default = '')
    %                 .AccessLevel (default = 1)
    % 
    % Output(s):
    %   this        - [obj] MEFSession_3p0 object
    %
    % See also get_sessinfo.

	% Copyright 2020 Richard J. Cui. Created: Thu 02/06/2020 10:07:26.965 AM
	% $Revision: 0.5 $  $Date: Thu 04/02/2020 10:41:14.141 AM $
	%
    % Multimodel Neuroimaging Lab (Dr. Dora Hermes)
    % Mayo Clinic St. Mary Campus
    % Rochester, MN 55905
    %
    % Email: richard.cui@utoronto.ca
    
    % =====================================================================
    % properties
    % =====================================================================
    % properties of importing session
    % -------------------------------
    properties
        MetaData            % session metadata (see read_mef_session_metadata_3p0.m
                            % for the detail)
        PathToSession       % not include session name and extension
        SessionName         % not include extension
        SessionExt          % session extension (includes the '.')
   end
    
    % =====================================================================
    % methods
    % =====================================================================
    % the constructor
    % ----------------
    methods
        function this = MEFSession_3p0(varargin)
            % MEFSession_3p0 Construct an instance of this class
            % ==================================================
            % parse inputs
            % -------------
            default_sp = ''; % default session path
            default_pw = struct('Level1Password', '', 'Level2Password', '',...
                'AccessLevel', 1);
            
            % parse rules
            p = inputParser;
            p.addOptional('filename', default_sp, @ischar);
            p.addOptional('password', default_pw, @isstruct)
            
            % parse and retrun the results
            p.parse(varargin{:});
            q = p.Results;
            filename = q.filename;
            password = q.password;
            
            % operations during construction
            % ------------------------------
            % initialize super classes
            this@MEFSession;
            this@MultiscaleElectrophysiologyFile_3p0;
            
            % set MEF version to serve
            if isempty(this.MEFVersion) == true
                this.MEFVersion = 3.0;
            elseif this.MEFVersion ~= 3.0
                error('MEFSession_3p0:invalidMEFVer',...
                    'invalid MEF version; this function can serve only MEF 3.0')
            end % if            
            
            % set session info
            [sesspath, channames] = this.findSessPath(filename);
            if ~isempty(sesspath)
                this.setSessionInfo(sesspath, password);
            end % if
            if ~isempty(channames)
                this.SelectedChannel = channames;
            end % if
        end % function
    end % methods
    
    % static methods
    % -------------
    methods (Static)
        
    end % methods

    % other methods
    % -------------
    methods
        [sesspath, channames] = findSessPath(this, filename) % find session path and channel name
        metadata = read_mef_session_metadata_3p0(this, varargin) % get session metadata of MEF 3.0
        valid_yn = checkSessValid(this, varargin) % check validity of session info
        [X, t] = import_sess(this, varargin) % import session of MEF 3.0 data
        metadata = setSessionInfo(this, sesspath, password) % set session information
    end % methods
end % classdef

% [EOF]