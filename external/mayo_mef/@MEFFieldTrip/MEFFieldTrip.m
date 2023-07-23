classdef MEFFieldTrip < handle
    % MEFFieldTrip process MEF in FieldTrip
    %  
    % Syntax:
    %
    % Input(s):
    %
    % Output(s):
    %
    % Note:
    %
    % See also .
    
    % Copyright 2020 Richard J. Cui. Created: Sat 03/21/2020 10:35:23.147 PM
    % $Revision: 0.2 $  $Date:  Wed 04/01/2020 11:55:41.288 PM $
    %
    % Multimodel Neuroimaging Lab (Dr. Dora Hermes)
    % Mayo Clinic St. Mary Campus
    % Rochester, MN 55905, USA
    %
    % Email: richard.cui@utoronto.ca


    % =====================================================================
    % properties
    % =====================================================================    
    properties
        FileType        % file type in field trip
    end % properties
    
    % =====================================================================
    % Constructor
    % =====================================================================
    methods
        function this = MEFFieldTrip()
            
        end % function
    end % methods
    
    % =====================================================================
    % methods
    % =====================================================================
    methods
        dc_event = findDiscontEvent(this, start_end, unit) % process discontinuity events
    end % methods
end % classdef

% [EOF]