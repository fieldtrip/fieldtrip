classdef MEFFieldTrip_3p0 < MEFSession_3p0 & MEFFieldTrip
    % MEFFIELDTRIP_3P0 process MEF 3.0 in FieldTrip
    % 
    % Syntax:
    %   this = MEFFieldTrip_3p0(filename)
    %   this = MEFFieldTrip_3p0(__, password)
    %
    % Input(s):
    %   filename    - [char] session path or channel path or dataset name
    %   password    - [struct] (opt) password (default: empty)
    %                 .Level1Password
    %                 .Level2Password
    %                 .AccessLevel
    %
    % Output(s):
    %   this        - [obj] MEFFieldTrip_3p0 object
    %
    % Note:
    %
    % See also .
    
    % Copyright 2020 Richard J. Cui. Created: Sat 03/21/2020 10:35:23.147 PM
    % $Revision: 0.2 $  $Date: Sun 03/22/2020  7:54:27.234 AM $
    %
    % Multimodel Neuroimaging Lab (Dr. Dora Hermes)
    % Mayo Clinic St. Mary Campus
    % Rochester, MN 55905
    %
    % Email: richard.cui@utoronto.ca
    
    % =====================================================================
    % properties
    % =====================================================================
    properties
        FileType        % file tpe in field trip
    end
    
    % =====================================================================
    % constructor
    % =====================================================================
    methods
        function this = MEFFieldTrip_3p0(varargin)
            % class constructor
            % =================
            % parse inputs
            % ------------
            % default
            default_sp = '';
            default_pw = struct('Level1Password', '',...
                'Level2Password', '', 'AccessLevel', 1);
            
            % parse rules
            p = inputParser;
            p.addOptional('sesspath', default_sp, @ischar);
            p.addOptional('password', default_pw, @isstruct);
            
            % parse the return the results
            p.parse(varargin{:});
            q = p.Results;
            sesspath = q.sesspath;
            password = q.password;
            
            % operations during construction
            % ------------------------------
            % call super class
            this@MEFFieldTrip;
            this@MEFSession_3p0;
            
            % set class properties
            this.FileType = 'mayo_mef30';
            
            % set session information
            if ~isempty(sesspath)
                this.setSessionInfo(sesspath, password);
            end % if
        end %function
    end % methods
    
    % =====================================================================
    % methods
    % =====================================================================
    methods (Static = true)

    end % 
    
    % other methods
    % -------------
    methods
        [sesspath, channames] = findSessPath(this, filename) % find session path and channel name
        hdr = getHeader(this, channames) % get header information of MEF 3.0 session
        evt = getEvent(this, channames) % get MEF 3.0 events for FieldTrip
        dat = getData(this, varargin) % read data from MEF 3.0 dataset for FieldTrip
    end % methods
end

% [EOF]