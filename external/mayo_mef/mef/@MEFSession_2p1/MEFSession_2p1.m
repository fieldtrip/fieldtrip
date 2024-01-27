classdef MEFSession_2p1 < MultiscaleElectrophysiologyFile_2p1 & MEFSession
	% Class MEFSESSION_2P1 processes MEF 2.1 session
    % 
    % Syntax:
    %   this = MEFSession_2p1
    %   this = MEFSession_2p1(filename)
    %   this = MEFSession_2p1(filename, password)
    %
    % Input(s):
    %   filename    - [str] (opt) MEF 2.1 session path or data file
    %   password    - [struct] (opt) structure of MEF 2.1 passowrd
    %                 .subject (default = '')
    %                 .session (default = '')
    %                 .data (default = '')
    % 
    % Output(s):
    %   this        - [obj] MEFSession_2p1 object
    %
    % See also get_sessinfo.

	% Copyright 2019-2020 Richard J. Cui. Created: Mon 12/30/2019 10:52:49.006 PM
	% $Revision: 1.4 $  $Date: Thu 04/02/2020 12:07:48.691 PM $
	%
	% 1026 Rocky Creek Dr NE
	% Rochester, MN 55906, USA
	%
	% Email: richard.cui@utoronto.ca

    % =====================================================================
    % properties
    % =====================================================================
    % properties of session information
    % ---------------------------------
    properties 

    end % properties
 
    % =====================================================================
    % methods
    % =====================================================================
    % the constructor
    % ----------------
    methods 
        function this = MEFSession_2p1(varargin)
            % MEFSession_2p1 contractor
            % =========================
            % parse inputs
            % ------------
            % defaults
            default_sp = ''; % default session path
            default_pw = struct('Subject', '', 'Session', '', 'Data', '');
            % parse rules
            p = inputParser;
            p.addOptional('filename', default_sp, @isstr);
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
            
            % set and check MEF version
            if isempty(this.MEFVersion) == true
                this.MEFVersion = 2.1;
            elseif this.MEFVersion ~= 2.1
                error('MEFSession_2p1:invalidMEFVer',...
                    'invalid MEF version; this function can serve only MEF 2.1')
            end % if            

            % set session info
            this.Password = password; % set password
            [sesspath, channames] = this.findSessPath(filename);
            if ~isempty(sesspath)
                this.SessionPath = sesspath; % set session path directory
                this.get_sessinfo;
            end % if
            if ~isempty(channames)
                this.SelectedChannel = channames;
            end % if
            
        end
    end % methods
    
    % static methods
    % -------------
    methods (Static)
        
    end % methods

    % other methods
    % -------------
    methods
        valid_yn = checkSessValid(this, varargin) % check validity of session info
        [sessionifo, unit] = get_info_data(this) % get session info from data
        [X, t] = import_sess(this, varargin) % import session of MEF 2.1 data
        [sesspath, channames] = findSessPath(this, filename) % find session path and channel name
    end % methods
    
end % classdef

% [EOF]
