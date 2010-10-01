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
      
          % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % missing data
      if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
     
      % multiple outputs
      if size(Y,2) > 1
        obj = ft_mv_noutput('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      obj.mu = nanmean(X);
      obj.sigma = nanstd(X);
      obj.sigma(obj.sigma==0) = 1; % bug fix
      
    end
    
    
    function Y = test(obj,X)
      
      Y = bsxfun(@minus,X,obj.mu);
      Y = bsxfun(@rdivide,Y,obj.sigma);
      
    end
    
    
  end
  
end
