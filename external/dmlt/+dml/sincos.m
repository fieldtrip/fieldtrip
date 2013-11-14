classdef sincos < dml.method
% SINCOS circular regression by decomposing into sine and cosine.
%
%   DESCRIPTION
%   Angle is represented as sine and cosine on which regularized linear
%   regression is performed
%
%   EXAMPLE
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

    properties        
        
      regressor = dml.enet('family','gaussian','alpha',1);
        
      sinreg; % regressor applied to sine
      cosreg; % regressor applied to cosine
      
    end
    
    methods
        
        function obj = sincos(varargin)
            
            obj = obj@dml.method(varargin{:});            
        end        
        
        function obj = train(obj,X,Y)
          
          obj.sinreg = obj.regressor.train(X,sin(Y));
          obj.cosreg = obj.regressor.train(X,cos(Y));
          
        end
        
        function Y = test(obj,X)
          
          Y = atan2(obj.sinreg.test(X),obj.cosreg.test(X));
          
        end
        
    end
    
end