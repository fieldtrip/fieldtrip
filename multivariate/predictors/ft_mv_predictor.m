classdef ft_mv_predictor < ft_mv_method
%FT_MV_PREDICTOR abstract class 
%
% Copyright (c) 2009, Marcel van Gerven
    
  methods

    function obj = ft_mv_predictor(varargin)

     obj@ft_mv_method(varargin{:});
      
    end
    
  end

end