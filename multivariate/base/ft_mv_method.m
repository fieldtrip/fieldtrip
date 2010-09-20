classdef ft_mv_method
% FT_MV_METHOD abstract handle class for multivariate methods
%
%   Copyright (c) 2009, Marcel van Gerven


  properties
  
    verbose = false;
    
  end
  
  methods
    
    function obj = ft_mv_method(varargin)
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        else
          error(sprintf('unrecognized fieldname %s',varargin{i}));
        end
      end
      
    end
    
    function [m,desc] = model(obj)
      % default behaviour when we ask for a model (override in subclass)
       
      m = {};
      desc = {};
      
    end
    
  end
  
  methods (Abstract)
      obj = train(obj,X,Y)
      Y = test(obj,X)
   end
  
end
