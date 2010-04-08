classdef whitener < preprocessor
%WHITENER whitens the data
%
% NOTE:
% data is assumed to be standardized!
%
% Copyright (c) 2009, Marcel van Gerven
%
   
    methods
      
      function obj = whitener(varargin)
        
        obj = obj@preprocessor(varargin{:});
        
      end
      
      function Y = map(obj,X)
        % whiten
        
        Y = X*obj.params.wmat';
          
      end
      
      function X = unmap(obj,Y)
        % unwhiten
        
        X = Y*obj.params.uwmat';
        
      end
      
      function p = estimate(obj,X,Y)
        
        [E, D] = eig(cov(X,1));
        
        p.wmat = sqrt(D) \ E';
        p.uwmat = E * sqrt(D);
        
      end
      
    end

end