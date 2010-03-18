classdef logreg < classifier
% logistic regression; various flavours using Mark Schmidt's L1General
% package
%
%   Copyright (c) 2010, Marcel van Gerven

  properties
    
    lambda = 0; % L1 regularization parameter
    nu = 0;     % L2 regularization parameter
    
  end

  methods
    
    function obj = logreg(varargin)
      
      obj = obj@classifier(varargin{:});
      
    end
    
    function p = estimate(obj,X,Y)
     
      X = [ones(size(X,1),1) X];
    
      nvars = size(X,2);
    
      opt.display = 0;
      
      nclasses = numel(unique(Y(:,1)));
      
      if nclasses == 2
        Y = 2*Y - 3; % convert to -1 / + 1
        funObj = @(w)LogisticLoss(w,X,Y);
      else
        funObj = @(w)SoftmaxLoss2(w,X,Y,nclasses);
      end
      
      if ~isempty(obj.params) && isfield(obj.params,'model') && ~isempty(obj.params.model)
        % initial model has been specified
        w_init = obj.params.model;
      else
        w_init = zeros(nvars,nclasses-1);
      end
      w_init = w_init(:);
      
      if obj.lambda == 0 && obj.nu == 0
        % unregularized logistic regression
        
        opt.Method = 'newton';
        p.model = minFunc(funObj,w_init,opt);
        
      elseif obj.lambda == 0 && obj.nu ~= 0
        % L2 regularized logistic regression
        
        nu = obj.nu*ones(nvars,nclasses-1);
        nu(1,:) = 0; % Don't regularize bias elements
        
        funObjL2 = @(w)penalizedL2(w,funObj,nu(:));
        
        opt.Method = 'newton';
        p.model = minFunc(funObjL2,w_init,opt);
        
      elseif obj.lambda ~= 0 && obj.nu == 0
        % L1 regularized logistic regression
        
        lambda = obj.lambda*ones(nvars,nclasses-1);
        lambda(1,:) = 0; % Don't regularize bias elements
        
        p.model = L1GeneralProjection(funObj,w_init,lambda(:));
        
      else
        % elastic net logistic regression
        
        nu = obj.nu*ones(nvars,nclasses-1);
        nu(1,:) = 0; % Don't regularize bias elements
        funObjL2 = @(w)penalizedL2(w,funObj,nu(:));
        
        % Set up regularizer
        lambda = obj.lambda*ones(nvars,nclasses-1);
        lambda(1,:) = 0; % Don't regularize bias elements
        
        p.model = L1GeneralProjection(funObjL2,w_init,lambda(:));
        
      end
      
      p.model = reshape(p.model,nvars,nclasses-1);
      
    end
    
    function Y = map(obj,X)
      
      Y = exp([ones(size(X,1),1) X] * [zeros(size(X,2)+1,1) obj.params.model]);
      Y = bsxfun(@rdivide,Y,sum(Y,2));
      
    end
    
    function [m,desc] = getmodel(obj)
      % return the parameters wrt a class label in some shape
      
      m =  mat2cell(obj.params.model(2:end,:)',ones(1,size(obj.params.model,2)),size(obj.params.model,1)-1);
      desc = {'unknown'};
      
    end
    
  end
end
