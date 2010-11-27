classdef ft_mv_rkls < ft_mv_kernelmethod
%FT_MV_RKLS regularized kernel least squares
%
% Copyright (c) 2008, Marcel van Gerven, Jason Farquhar
  
  methods
    
    function obj = ft_mv_rkls(varargin)
      
      obj = obj@ft_mv_kernelmethod(varargin{:});
      
    end
    
    function [weights,f,J] = estimate(obj,K,Y)
      % this kernelmethod's estimation function
      
       if ~isempty(obj.weights)
        [weights,f,J] = rkls(K,Y,obj.C,'verb',-1,'alphab',obj.weights);
      else
        [weights,f,J] = rkls(K,Y,obj.C,'verb',-1);
      end
      
    end
    
  end
end
