classdef ft_mv_timeseries < ft_mv_method
%FT_MV_TIMESERIES abstract class 
%
% Copyright (c) 2009, Marcel van Gerven
    
  properties

   indims % input dimensions; required to retrieve the trials x features x time structure
   
  end

  methods

    function obj = ft_mv_timeseries(varargin)

     obj@ft_mv_method(varargin{:});
      
    end
    
    function Z = predict(obj,X)
      % return predictions instead of posteriors for discrete variables
      
      Z = obj.test(X);
      
      if all(sum(Z,2) == 1)
         % convert posteriors into classifications
        [tmp,Z] = max(Z,[],2);
      end
      
    end
    
  end

end