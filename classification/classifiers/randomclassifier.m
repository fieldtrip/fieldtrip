classdef randomclassifier < classifier
%RANDOMCLASSIFIER returns a random classification
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: randomclassifier.m,v $
%

  properties
    nclasses;
  end
  
  methods
    function obj = randomclassifier(varargin)
      obj = obj@classifier(varargin{:});
    end
    
    function obj = train(obj,data,design)
      % does nothing
      
      obj.nclasses = design.nunique;
      
    end
    
    function post = test(obj,data)
      
      data = data.collapse();
      
      % random classification
      post = rand(size(data,1),obj.nclasses);
      post = post ./ repmat(sum(post,2),[1 size(post,2)]);
      
      post = dataset(post);
      
    end
    
  end
end
