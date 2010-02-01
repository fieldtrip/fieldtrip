classdef standardizer < preprocessor
%STANDARDIZER standardizes data to have mean 0 and standard deviation 1
%
% SEE ALSO
%   zscore.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: standardizer.m,v $
%

    methods
        
      function obj = standardizer(varargin)
        
        obj = obj@preprocessor(varargin{:});
        
      end
                 
      function M = map(obj,U)
        
        M = bsxfun(@minus,U.X,obj.params.mu);
        M = dataset(bsxfun(@rdivide,M,obj.params.sigma));
        
      end
      
      function U = unmap(obj,M)
        
         U = bsxfun(@times,M.X,obj.params.sigma);
         U = dataset(bsxfun(@plus,U,obj.params.mu));
       
      end
     
      function p = estimate(obj,data,design)
        
        p.mu = mynanmean(data.X);
        p.sigma = mynanstd(data.X);
        p.sigma(p.sigma==0) = 1; % bug fix
        
      end
      
    end
    
end
