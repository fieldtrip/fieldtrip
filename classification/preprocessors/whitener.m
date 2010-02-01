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
      
      function M = map(obj,U)
        % whiten
        
        M = dataset(U.X*obj.params.wmat');
          
      end
      
      function U = unmap(obj,M)
        % unwhiten
        
        U = dataset(M.X*obj.params.uwmat');
        
      end
      
      function p = estimate(obj,data,design)
        
        [E, D] = eig(cov(data.X,1));
        
        p.wmat = sqrt(D) \ E';
        p.uwmat = E * sqrt(D);
        
      end
      
    end

end