classdef ft_mv_glmlogreg < ft_mv_predictor
% logistic regression; various flavours using Friedman's glmnet package
% data is assumed to be standardized; so no bias term is needed
%
%   Copyright (c) 2010, Marcel van Gerven

  properties
    
    lambda = 0.1; % L1 parameter
    alpha = 1; % mixing parameter; alpha < 1 = elastic net
    
    weights   % regression coefficients
  end

  methods
    
    function obj = ft_mv_glmlogreg(varargin)
      
      obj = obj@ft_mv_predictor(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
     
      nclasses = numel(unique(Y(:,1)));
      
      opts = glmnetSet;
      
      opts.alpha = obj.alpha; % mixing parameter
      opts.lambda = obj.lambda; 

      if nclasses == 2
  
        fit=glmnet(X,Y,'binomial',opts);
         
      else
      
        fit=glmnet(X,Y,'multinomial',opts);
        % not yet implemented
        
      end
      
      obj.weights = fit.beta;
      
    end
    
    function Y = test(obj,X)
      
      Y = exp(X * [obj.weights zeros(size(X,2),1)]);
      Y = 1 - bsxfun(@rdivide,Y,sum(Y,2));
      
    end
    
    function [m,desc] = model(obj)
      % return the parameters wrt a class label in some shape
      
      % weight-vector for class 1 is always zero
      m =  mat2cell(obj.weights',ones(1,size(obj.weights,2)),size(obj.weights,1));
      m{length(m)+1,1} = zeros(size(m{1}));
      
      desc = cell(length(m),1);
      for j=1:length(m)
        desc{j} = sprintf('regression coefficients for class %d',j);
      end
      
    end
    
  end
end
