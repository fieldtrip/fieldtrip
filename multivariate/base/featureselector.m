classdef featureselector < mvmethod
%FEATURESELECTOR featureselector method class
%
% During operation, the featureselector takes data and
% produces a reduced dataset with (typically) a smaller number of features.
%
% PARAMETERS:
%   'subset'    : the subset which can be set manually using this class
%
% SEE ALSO
% doc featureselectors
%
%   Copyright (c) 2008, Marcel van Gerven
   
  properties
    
    procedure             % the multivariate analysis used for prediction
      
    validator             % the used validator
      
    metric = 'accuracy';  % evaluation metric
  
  end
  

  methods
    
    function obj = featureselector(varargin)
      
      % parse options
      for i=1:2:length(varargin)
        if ismember(varargin{i},fieldnames(obj))
          obj.(varargin{i}) = varargin{i+1};
        end
      end
            
    end
    
    function Y = map(obj,X)
      
      if obj.verbose
        fprintf('using %d out of %d features\n',numel(obj.params.subset),size(X,2));
      end
      
      Y = X(:,obj.params.subset);
      
    end
    
    function [model,desc] = getmodel(obj)
      % return used subset
      
      model = {obj.params.subset};
      desc = {'indices of the used features'};
      
    end
    
  end
  
end
