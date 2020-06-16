classdef MEFFieldTrip_2p1 < MEFSession_2p1 & MEFFieldTrip
    % MEFFIELDTRIP_3P0 process MEF 2.1 in FieldTrip
    % 
    % Syntax:
    %   this = MEFFieldTrip_2p1(filename)
    %   this = MEFFieldTrip_2p1(__, password)
    %
    % Input(s):
    %   filename    - [char] session path or dataset name
    %   password    - [struct] (opt) structure of MEF 2.1 passowrd
    %                 .subject (default = '')
    %                 .session (default = '')
    %                 .data (default = '')
    %
    % Output(s):
    %   this        - [obj] MEFFieldTrip_2p1 object
    %
    % Note:
    %
    % See also .
    
    % Copyright 2020 Richard J. Cui. Created: Wed 04/01/2020 11:55:41.288 PM
    % $Revision: 0.1 $  $Date: Wed 04/01/2020 11:55:41.288 PM $
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
        function this = MEFFieldTrip_2p1(varargin)
            % class constructor
            % =================
            % parse inputs
            % ------------
            
            % operations during construction
            % ------------------------------
            % call super class
            this@MEFFieldTrip;
            this@MEFSession_2p1(varargin{:});
            
            % set class properties
            this.FileType = 'mayo_mef21';
            
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
        hdr = getHeader(this, channames) % get header information of MEF 2.1 session
        evt = getEvent(this, channames) % get MEF 2.1 events for FieldTrip
        dat = getData(this, varargin) % read data from MEF 2.1 dataset for FieldTrip
    end % methods
end

% [EOF]