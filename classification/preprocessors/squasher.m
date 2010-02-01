classdef squasher < preprocessor
% SQUASHER squashes data between 0 and 1
%
% PARAMETERS:
%  dmin; % minimum of training data
%  dmax; % maximum of training data
%    
%   Copyright (c) 2009, Marcel van Gerven

  methods
    
    function obj = squasher(varargin)
      
      obj = obj@preprocessor(varargin{:});
      
    end
    
    function M = map(obj,U)
      
      M = bsxfun(@minus,U.X, obj.params.dmin);
      M = dataset(bsxfun(@rdivide,M, obj.params.dmax));
      
    end
    
    function U = unmap(obj,M)
      
      U = bsxfun(@times,M.X, obj.params.dmax);
      U = dataset(bsxfun(@plus,U, obj.params.dmin));
      
    end
    
    function p = estimate(obj,data,design)
      
      p.dmin = min(data.X);
      X = bsxfun(@minus,data.X, obj.dmin);
      p.dmax = max(X);
      p.dmax(p.dmax == 0) = 1;
      
    end
    
  end
end
