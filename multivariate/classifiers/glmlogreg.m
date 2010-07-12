classdef glmlogreg < classifier
% logistic regression; various flavours using Friedman's glmnet package
% data is assumed to be standardized; so no bias term is needed
%
%   Copyright (c) 2010, Marcel van Gerven

  properties
    
    lambda = 0.1; % L1 parameter
    alpha = 1; % mixing parameter; alpha < 1 = elastic net
    
  end

  methods
    
    function obj = glmlogreg(varargin)
      
      obj = obj@classifier(varargin{:});
      
    end
    
    function p = estimate(obj,X,Y)
     
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
      
      p.model = fit.beta;
      
    end
    
    function Y = map(obj,X)
      
      Y = exp(X * [obj.params.model zeros(size(X,2),1)]);
      Y = 1 - bsxfun(@rdivide,Y,sum(Y,2));
      
    end
    
    function [m,desc] = getmodel(obj)
      % return the parameters wrt a class label in some shape
      
      % weight-vector for class 1 is always zero
      m =  mat2cell(obj.params.model',ones(1,size(obj.params.model,2)),size(obj.params.model,1));
      m{length(m)+1,1} = zeros(size(m{1}));
      
      desc = cell(length(m),1);
      for j=1:length(m)
        desc{j} = sprintf('regression coefficients for class %d',j);
      end
      
    end
    
  end
end
