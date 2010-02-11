classdef standardizer < preprocessor
%STANDARDIZER standardizes data to have mean 0 and standard deviation 1
%
% SEE ALSO
%   zscore.m
%
%   Copyright (c) 2008, Marcel van Gerven


    methods
        
      function obj = standardizer(varargin)
        
        obj = obj@preprocessor(varargin{:});
        
      end
                 
      function M = map(obj,U)
        
        M = bsxfun(@minus,U,obj.params.mu);
        M = bsxfun(@rdivide,M,obj.params.sigma);
        
      end
      
      function U = unmap(obj,M)
        
         U = bsxfun(@times,M,obj.params.sigma);
         U = bsxfun(@plus,U,obj.params.mu);
       
      end
     
      function p = estimate(obj,data,design)
        
        p.mu = mynanmean(data);
        p.sigma = mynanstd(data);
        p.sigma(p.sigma==0) = 1; % bug fix
        
      end
      
    end
    
end
