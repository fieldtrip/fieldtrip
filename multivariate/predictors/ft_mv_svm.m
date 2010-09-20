classdef ft_mv_svm < ft_mv_kernelmethod
%FT_MV_SVM support vector machine
%
% Copyright (c) 2008, Marcel van Gerven, Jason Farquhar
  
  methods
    
    function obj = ft_mv_svm(varargin)
      
      obj = obj@ft_mv_kernelmethod(varargin{:});
      
    end
    
    function [weights,f,J] = estimate(obj,K,Y)
      % this kernelmethod's estimation function
      
      [weights,f,J] = l2svm_cg(K,Y,obj.C,'verb',-1);
      
    end
    
  end
end
