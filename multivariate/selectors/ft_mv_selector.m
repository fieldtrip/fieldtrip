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
