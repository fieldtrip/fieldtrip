classdef ft_mv_klr < ft_mv_kernelmethod
%FT_MV_KLR kernel logistic regression
%
% Copyright (c) 2008, Marcel van Gerven, Jason Farquhar
  
  methods
    
    function obj = ft_mv_klr(varargin)
      
      obj = obj@ft_mv_kernelmethod(varargin{:});
      
    end
    
    function [weights,f,J] = estimate(obj,K,Y)
      % this kernelmethod's estimation function
      
      [weights,f,J] = klr_cg(K,Y,obj.C,'verb',-1);
      
    end
    
  end
end
