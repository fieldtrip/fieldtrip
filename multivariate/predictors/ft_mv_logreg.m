classdef ft_mv_logreg < ft_mv_predictor
% FT_MV_LOGREG logistic regression; various flavours using Mark Schmidt's L1General
% package; elastic net implementation based on Friedman et al
%
%   Copyright (c) 2010, Marcel van Gerven, Ali Bahramisharif

  properties
    
    L1 = 0; % L1 regularization parameter
    L2 = 0; % L2 regularization parameter
    
    tol; % convergence criterion of elastic net
    
    weights; % regression coefficients
    
  end

  methods
    
    function obj = ft_mv_logreg(varargin)
      
      obj = obj@ft_mv_predictor(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)
     
      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % missing data
      if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
     
      % multiple outputs
      if size(Y,2) > 1
        obj = ft_mv_noutput('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      X = [ones(size(X,1),1) X];
    
      nvars = size(X,2);
    
      opt.display = 0;
      
      nclasses = max(Y(:));
      
      if nclasses == 2
        % convert to -1 / + 1
        funObj = @(w)LogisticLoss(w,X,3 - 2*Y);
      else
        funObj = @(w)SoftmaxLoss2(w,X,Y,nclasses);
      end
      
      if ~isempty(obj.weights)
        % initial model has been specified
        w_init = obj.weights;
      else

        % class 1 only
        if nclasses == 1
          obj.weights = zeros(nvars,1);
          obj.weights(1) = 1;
          return;
        end
        
        % deal with only one class as input
        u = unique(Y);
        if length(u) == 1
          obj.weights = zeros(nvars,nclasses-1);
          obj.weights(1,u-1) = -1;
          return;
        end
        
        w_init = zeros(nvars,nclasses-1);

      end
      w_init = w_init(:);
      
      if obj.L1 == 0 && obj.L2 == 0
        % unregularized logistic regression
        
        opt.Method = 'newton';
        obj.weights = minFunc(funObj,w_init,opt);
        
      elseif obj.L1 == 0 && obj.L2 ~= 0
        % L2 regularized logistic regression
        
        L2 = obj.L2*ones(nvars,nclasses-1);
        L2(1,:) = 0; % Don't regularize bias elements
        
        funObjL2 = @(w)penalizedL2(w,funObj,L2(:));        
        opt.Method = 'newton';
        obj.weights = minFunc(funObjL2,w_init,opt);
        
      elseif obj.L1 ~= 0 && obj.L2 == 0
        % L1 regularized logistic regression
        
        L1 = obj.L1*ones(nvars,nclasses-1);
        L1(1,:) = 0; % Don't regularize bias elements
        
        obj.weights = L1GeneralProjection(funObj,w_init,L1(:));
        
      else
        % elastic net logistic regression
        
        % Mark Schmidt implementation
        %         L2 = obj.L2*ones(nvars,nclasses-1);
        %         L2(1,:) = 0; % Don't regularize bias elements
        %         funObjL2 = @(w)penalizedL2(w,funObj,L2(:));
        %
        %         % Set up regularizer
        %         L1 = obj.L1*ones(nvars,nclasses-1);
        %         L1(1,:) = 0; % Don't regularize bias elements
        %
        %         p.model = L1GeneralProjection(funObjL2,w_init,L1(:));

        opt.L1=obj.L1;            
        opt.L2=obj.L2;                 
        opt.verbose=obj.verbose; 
        opt.tol=obj.tol;
        
        obj.weights = -elasticlr(X,Y,opt);
         
      end
      
      obj.weights = reshape(obj.weights,nvars,nclasses-1);
      
    end
    
    function Y = test(obj,X)
      
      Y = exp([ones(size(X,1),1) X] * [obj.weights zeros(size(X,2)+1,1)]);
      Y = bsxfun(@rdivide,Y,sum(Y,2));
      
    end
    
    function [m,desc] = model(obj)
      % return the parameters wrt a class label in some shape
      
      % weight-vector for class 1 is always zero
      m =  mat2cell(obj.weights(2:end,:)',ones(1,size(obj.weights,2)),size(obj.weights,1)-1);
      m{length(m)+1,1} = zeros(size(m{1}));
      
      desc = cell(length(m),1);
      for j=1:length(m)
        desc{j} = sprintf('regression coefficients for class %d',j);
      end
      
    end
    
  end
end
