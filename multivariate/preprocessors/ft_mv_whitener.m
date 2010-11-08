classdef ft_mv_whitener < ft_mv_standardizer
%WHITENER standardizes and whitens the data
%
% Copyright (c) 2009, Marcel van Gerven
   
  properties
    
    wmat  % whitening matrix
    
  end

  methods
    
    function obj = ft_mv_whitener(varargin)
      
      obj = obj@ft_mv_standardizer(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
      
      if nargin<3, Y = []; end
      
      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % missing data
      if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
     
      obj = train@ft_mv_standardizer(obj,X,Y);
      
      if obj.verbose, fprintf('whitening data\n'); end
      
      [E, D] = eig(cov(X,1));
      obj.wmat = sqrt(D) \ E';
      
    end
    
    function Y = test(obj,X)
      % whiten
      
      Y = test@ft_mv_standardizer(obj,X);
      
      Y = Y*obj.wmat';
      
    end
    
    function X = invert(obj,Y)
      % invert mapping
      
      X = Y*inv(obj.wmat)';
      
      X = invert@ft_mv_standardizer(obj,X);
      
    end
    
  end
  
end