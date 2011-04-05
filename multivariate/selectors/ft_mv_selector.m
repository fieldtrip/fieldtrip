classdef ft_mv_selector < ft_mv_method
%FT_MV_SELECTOR featureselector method class
%
% The featureselector takes data and produces a reduced dataset with a smaller number of features.
%
%   Copyright (c) 2008, Marcel van Gerven
   
  properties
    
    validator % the validator used to determine the final used feature subset
        
    subset % the used feature subset
    
  end
  

  methods
    
    function obj = ft_mv_selector(varargin)
      
      obj = obj@ft_mv_method(varargin{:});
      
%       if isempty(obj.validator)
%         fprintf('validator not specified, using 5-fold CV with SVM as classifier and accuracy as performance metric\n');
%         obj.validator = ft_mv_crossvalidator('nfolds',5,'mva',ft_mv_svm,'metric','accuracy');
%       end
                  
    end
    
    function obj = train(obj,X,Y)
      % in case the subset is manually specified
      
    end
    
    function Y = test(obj,X)
            
      Y = X(:,obj.subset);
      
    end
    
    function [model,desc] = model(obj)
      % return logical array with ones indicating selected features
      
      model = zeros(obj.indims);
      model(obj.params.subset) = 1;      
      model = {model};
      desc = {'logical array with selected features'};
      
    end
    
  end
  
end
