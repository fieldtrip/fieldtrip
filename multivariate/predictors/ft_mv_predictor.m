classdef ft_mv_predictor < ft_mv_method
%FT_MV_PREDICTOR abstract class 
%
% Copyright (c) 2009, Marcel van Gerven
    
  methods

    function obj = ft_mv_predictor(varargin)

     obj@ft_mv_method(varargin{:});
      
    end
    
    function Z = predict(obj,X)
      % return predictions (maximum of the returned outputs)
      
      Z = obj.test(X);
      
      if size(Z,2) > 1
         % convert posteriors into classifications
        [tmp,Z] = max(Z,[],2);
      end
      
    end
    
  end

end