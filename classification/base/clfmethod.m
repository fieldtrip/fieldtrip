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
          
          if iscell(data)
            dims = cell(1,length(data));
            for c=1:length(data)
              dims{c} = data{c}.dims;
            end
          else
            dims = data.dims;
          end
          
          data = obj.map(data);
          
           if iscell(data)
            dims = cell(1,length(data));
            for c=1:length(data)
              data{c}.dims = dims{c};
              data{c}.ndims = length(dims{c});
            end
          else
            data.dims = dims;
            data.ndims = length(dims);
          end
          
        end
      end
      
      function data = untest(obj,data)
        % invert the mapping
        
        dims = data.dims;
        
        if iscell(data) && ~isa(obj,'transfer_learner')
          
          params = obj.params;
          
          for c=1:length(data)
            
            obj.params = params{c};
            data{c} = obj.test(data{c});
          end
          
        else
          
         if iscell(data)
            dims = cell(1,length(data));
            for c=1:length(data)
              dims{c} = data{c}.dims;
            end
          else
            dims = data.dims;
          end
          
          data = obj.unmap(data);
          
           if iscell(data)
            dims = cell(1,length(data));
            for c=1:length(data)
              data{c}.dims = dims{c};
              data{c}.ndims = length(dims{c});
            end
          else
            data.dims = dims;
            data.ndims = length(dims);
          end
          
          
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