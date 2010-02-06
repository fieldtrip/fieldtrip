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
      
      params; % the used parameters for the mapping/unmapping
 
    end
    
    methods                  
      
      function [m,desc] = getmodel(obj)
        % default behaviour when we ask for a model (override in subclass)
        
        if obj.verbose
          fprintf('don`t know how to return model for object of type %s; returning empty model and description\n',class(obj));
        end
        
        m = {};
        desc = {};
        
      end
      
      function obj = train(obj,data,design)
        
        if iscell(data) && ~isa(obj,'transfer_learner')
          
          params = cell(1,length(data));
          
          for c=1:length(data)
            obj = obj.train(data{c},design{c});
            params{c} = obj.params;
          end
          
          obj.params = params;
          
        else
          
          if iscell(data) && isa(obj,'transfer_learner') && ...
              length(unique(cellfun(@(x)(x.nfeatures),data))) > 1
            % check if datasets have the same number of features
              error('datasets must have the same number of features for transfer learning');
          end
          
          obj.params = obj.estimate(data,design);
          
        end
      end
      
      function data = test(obj,data)
        
        if iscell(data) && ~isa(obj,'transfer_learner')
          
          params = obj.params;
          
          for c=1:length(data)
            
            obj.params = params{c};
            data{c} = obj.test(data{c});
          end
                    
        else
          
          data = obj.map(data);
          
        end
      end
      
      function data = untest(obj,data)
        % invert the mapping
        
        if iscell(data) && ~isa(obj,'transfer_learner')
          
          params = obj.params;
          
          for c=1:length(data)
            
            obj.params = params{c};
            data{c} = obj.test(data{c});
          end
          
        else
          
          data = obj.unmap(data);
          
        end
      end
      
      function U = unmap(obj,M)
        % sometimes the inverse mapping does not exist
        
        error('inverse mapping does not exist for class %s',class(obj));
        
      end
      
    end
    
    methods(Abstract)
      
      M = map(obj,U)   % mapping function
      p = estimate(D); % parameter estimation
    end
    
end