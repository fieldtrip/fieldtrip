classdef MEFFieldTrip_3p0 < MEFSession_3p0 & MEFFieldTrip
    % MEFFIELDTRIP_3P0 process MEF 3.0 in FieldTrip
    % 
    % Syntax:
    %   this = MEFFieldTrip_3p0
    %   this = MEFFieldTrip_3p0(filename)
    %   this = MEFFieldTrip_3p0(__, password)
    %   this = MEFFieldTrip_3p0(__, 'SortChannel', sortchannel)
    %
    % Input(s):
    %   filename    - [char] (opt) session path or channel path or dataset name
    %   password    - [struct] (opt) password (default: empty)
    %                 .Level1Password
    %                 .Level2Password
    %                 .AccessLevel
    %   sortchannel - [char] (para) sort channel according to either 'alphabet' of
    %                 the channel names or 'number' of the acquisiton
    %                 channel number (default = 'alphabet')
    %
    % Output(s):
    %   this        - [obj] MEFFieldTrip_3p0 object
    %
    % Note:
    %
    % See also .
    
    % Copyright 2020 Richard J. Cui. Created: Sat 03/21/2020 10:35:23.147 PM
    % $Revision: 0.4 $  $Date: Thu 04/02/2020 12:27:12.159 AM $
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
            
            % operations during construction
            % ------------------------------
            % call super class
            this@MEFFieldTrip;
            this@MEFSession_3p0(varargin{:});
            
            % set class properties
            this.FileType = 'mayo_mef30';
            
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
        hdr = getHeader(this, channames) % get header information of MEF 3.0 session
        evt = getEvent(this, channames) % get MEF 3.0 events for FieldTrip
        dat = getData(this, varargin) % read data from MEF 3.0 dataset for FieldTrip
    end % methods
end

% [EOF]