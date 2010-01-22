classdef clfmethod
% clfmethod base class for classification toolbox methods
%   
%   This base class contains common properties
%   which may be called by all child methods
%
%   mainly deals with data handling
%
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: clfmethod.m,v $
%

    properties
      
      verbose = false;
      
    end
    
    methods                  
      
      function m = getmodel(obj)
        % default behaviour when we ask for a model (override in subclass)
        
        if obj.verbose
          fprintf('don`t know how to return model for object of type %s; returning empty model\n',class(obj));
        end
        
        m = [];
        
      end
      
      function obj = train(obj,data,design)
      % do nothing by default
      
      end
      
      function post = test(obj,data)
        % do nothing
        post = data;
      end
   
    end    
      
end