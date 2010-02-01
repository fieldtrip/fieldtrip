classdef randomclassifier < classifier
%RANDOMCLASSIFIER returns a random classification
%
%   Copyright (c) 2008, Marcel van Gerven

  
  methods
    
    function obj = randomclassifier(varargin)
      obj = obj@classifier(varargin{:});
    end
    
    function p = estimate(obj,data,design)
      % does nothing
      
      p.nclasses = design.nunique;
      
    end
    
    function post = map(obj,data)
      
      % random classification
      post = rand(data.nsamples,obj.params.nclasses);
      post = post ./ repmat(sum(post,2),[1 size(post,2)]);
      
      post = dataset(post);
      
    end
    
  end
end
