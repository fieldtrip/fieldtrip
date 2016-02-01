classdef standardizer < dml.method
% STANDARDIZER takes zscores.
% 
%   DESCRIPTION
%   Takes zscores such that data has mean 0 and standard deviation 1
%  
%   EXAMPLE
%   X = rand(10,3);
%   m = dml.standardizer;
%   m = m.train(X);
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
    mu    % mean
    sigma % standard deviation

  end
  
  methods
    
    function obj = standardizer(varargin)
      
      obj = obj@dml.method(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
      
      % handle multiple datasets
      if iscell(X)
        obj = dml.ndata('method',obj);
        obj = obj.train(X);
        return;
      end
      
      if obj.verbose, fprintf('standardizing data\n'); end
           
      obj.mu = nanmean(X);
      obj.sigma = nanstd(X);
            
    end
    
    
    function Y = test(obj,X)
      
      Y = X;
      
      idx = ~isnan(obj.mu);
      if any(idx)
        Y(:,idx) = bsxfun(@minus,X(:,idx),obj.mu(idx));
      end
      
      % columns with mean of nan are set to a small random value
      if any(~idx)
        Y(:,~idx) = 1e-6*randn(size(Y(:,~idx))); 
      end

      idx = ~isnan(obj.sigma) & ~(obj.sigma==0);
      if any(idx)
        Y(:,idx) = bsxfun(@rdivide,Y(:,idx),obj.sigma(idx));
      end
      
      % columns with mean of nan are set to a small random value
      if any(~idx)
        Y(:,~idx) = 1e-6*randn(size(Y(:,~idx))); 
      end
      
    end
    
    function X = invert(obj,Y)
    % INVERT inverts the standardization
    
      X = bsxfun(@times,Y,obj.sigma);
      X = bsxfun(@plus,X,obj.mu);
   
      if any(isnan(X(:))), warning('NaNs encountered'); end
      
    end
    
    
  end
  
end
