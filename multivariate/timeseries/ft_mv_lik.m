classdef ft_mv_lik < ft_mv_timeseries
%FT_MV_LIK classification by comparing mva likelihoods
% 
% We assume here that X(j,f,t) contains sequence j for feature f at time t
% Hence, all trials are of equal length. Y indicates the condition to which
% the sequence belongs
%
% Copyright (c) 2010, Marcel van Gerven
  
  properties

   mva % used mva
    
   indims % input dimensions; required when CV removes dimensions
   
  end

  methods
    
    function obj = ft_mv_lik(varargin)
      
      obj = obj@ft_mv_timeseries(varargin{:});
      
      
    end

    function obj = train(obj,X,Y)

      if ~isempty(obj.indims)
        X = reshape(X,[size(X,1) obj.indims]);
      end
      
      nclasses = max(Y(:));
      
      obj.mva = repmat({obj.mva},[1 nclasses]);
      for j=1:nclasses
        if obj.verbose, fprintf('estimating mva %d of %d\n',j,nclasses); end
        obj.mva{j} = obj.mva{j}.train(X(Y==j,:,:));
      end
      
      
    end
        
    function post = test(obj,X)   
      
      if ~isempty(obj.indims)
        X = reshape(X,[size(X,1) obj.indims]);
      end
      
      % get mva likelihoods
      for c=1:length(obj.mva)
        post(:,c) = obj.mva{c}.likelihood(X);
      end
      
    end
    
    function [m,d] = model(obj)
      
      m = [];
      d = [];
      
    end
  
  end
  
end
