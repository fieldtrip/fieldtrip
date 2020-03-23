classdef MEFSession_3p0 < MEFSession & MultiscaleElectrophysiologyFile_3p0
    % Class MEFSESSION_3P0 process MEF 3.0 session
    % 
    % Syntax:
    %   this = MEFSession_3p0
    %   this = __(sesspath)
    %   this = __(sesspath, password)
    %
    % Input(s):
    %   sesspath    - [str] (opt) MEF 3.0 session path (default = '')
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
	% $Revision: 0.4 $  $Date: Sat 03/21/2020 10:35:23.147 PM $
	%
	% 1026 Rocky Creek Dr NE
	% Rochester, MN 55906, USA
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
            p.addOptional('sesspath', default_sp, @isstr);
            p.addOptional('password', default_pw, @isstruct)
            
            % parse and retrun the results
            p.parse(varargin{:});
            q = p.Results;
            sesspath = q.sesspath;
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
            if ~isempty(sesspath)
                this.setSessionInfo(sesspath, password);
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
        metadata = read_mef_session_metadata_3p0(this, varargin) % get session metadata of MEF 3.0
        valid_yn = checkSessValid(this, varargin) % check validity of session info
        [X, t] = import_sess(this, varargin) % import session of MEF 3.0 data
        setSessionInfo(this, sesspath, password) % set session information
    end % methods
end % classdef

% [EOF]