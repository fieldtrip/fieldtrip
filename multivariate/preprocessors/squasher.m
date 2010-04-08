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
    
    function Y = map(obj,X)
      
      Y = bsxfun(@minus,X, obj.params.dmin);
      Y = bsxfun(@rdivide,Y, obj.params.dmax);
      
    end
    
    function X = unmap(obj,Y)
      
      X = bsxfun(@times,Y, obj.params.dmax);
      X = bsxfun(@plus,X, obj.params.dmin);
      
    end
    
    function p = estimate(obj,X,Y)
      
      p.dmin = min(X);
      X = bsxfun(@minus,X, p.dmin);
      p.dmax = max(X);
      p.dmax(p.dmax == 0) = 1;
      
    end
    
  end
end
