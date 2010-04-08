classdef randomclassifier < classifier
%RANDOMCLASSIFIER returns a random classification
%
%   Copyright (c) 2008, Marcel van Gerven

  
  methods
    
    function obj = randomclassifier(varargin)
      obj = obj@classifier(varargin{:});
    end
    
    function p = estimate(obj,X,Y)
      % does nothing
      
      p.nclasses = obj.nunique(Y);
      
    end
    
    function Y = map(obj,X)
      
      % random classification
      Y = rand(size(X,1),obj.params.nclasses);
      Y = bsxfun(@rdivide,Y,sum(Y,2));
     
    end
    
  end
end
