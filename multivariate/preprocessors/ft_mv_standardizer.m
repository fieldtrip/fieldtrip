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
      
      if nargin<3, Y = []; end;
      
      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % missing data
      if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
           
      obj.mu = nanmean(X);
      obj.sigma = nanstd(X);
      obj.sigma(obj.sigma==0) = 1; % bug fix
      
    end
    
    
    function Y = test(obj,X)
      
      Y = bsxfun(@minus,X,obj.mu);
      Y = bsxfun(@rdivide,Y,obj.sigma);
      
    end
    
    function X = invert(obj,Y)

      X = bsxfun(@times,Y,obj.sigma);
      X = bsxfun(@plus,X,obj.mu);
      
    end
    
    
  end
  
end
