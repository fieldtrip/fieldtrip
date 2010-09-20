classdef ft_mv_standardizer < ft_mv_preprocessor
%FT_MV_STANDARDIZER standardizes data to have mean 0 and standard deviation 1
%
% SEE ALSO
%   zscore.m
%
%   Copyright (c) 2008, Marcel van Gerven

  properties
    
    mu    % mean
    sigma % standard deviation

  end
  
  methods
    
    function obj = ft_mv_standardizer(varargin)
      
      obj = obj@ft_mv_preprocessor(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
      
      if iscell(X)
        obj.mu = cell(1,length(X));
        obj.sigma = cell(1,length(X));        
        for c=1:length(X)
        obj.mu{c} = nanmean(X{c});
        obj.sigma{c} = nanstd(X{c});
        obj.sigma{c}(obj.sigma{c}==0) = 1; % bug fix
        end
      else
        obj.mu = nanmean(X);
        obj.sigma = nanstd(X);
        obj.sigma(obj.sigma==0) = 1; % bug fix
      end
      
    end
    
    
    function Y = test(obj,X)
      
      if iscell(X)
        for c=1:length(X)
          Y{c} = bsxfun(@minus,X{c},obj.mu{c});
          Y{c} = bsxfun(@rdivide,Y{c},obj.sigma{c});
        end
      else
        Y = bsxfun(@minus,X,obj.mu);
        Y = bsxfun(@rdivide,Y,obj.sigma);
      end
      
    end
    
    
  end
  
end
