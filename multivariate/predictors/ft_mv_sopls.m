classdef ft_mv_sopls < ft_mv_predictor
% FT_MV_SOPLS sparse orthonormalized partial least squares
%
% lambda_min specifies the fraction of the maximal lambda value (computed
% by glmnet) for which the output will be returned. 
%
% alpha specifies the balance between ridge regression (alpha=0) and lasso
% regression (alpha=1).
%
% The number of steps taken to traverse the regularization path is given by nlambda.
% 
% if nlambda=0 then standard partial least squares will be used.
%
% pmax specifies the maximum number of variables included in the solution.
% 
% nhidden specifies the number of used components. If it is decreased
% during testing then the restricted number of components will be used
% (only in case of sparse partial least squares!)
% 
% load 69digits;
% p = ft_mv_sopls('nhidden',3,'verbose',true);
% p = p.train(X(1:100,:),I(1:100,:));
% r = p.test(X(101:111,:));
% images(r,1:10,[2 5]);
% figure
% images(I,101:110,[2 5]);
%
% Copyright (c) 2010, Marcel van Gerven, Tom Heskes


  properties
  
    nhidden = 5;      % number of hidden units; estimated using leave-one-out if an array
    
    A             % output matrix
    B             % input matrix
    C             % hidden bias vector
    
    score         % leave-one-out score
    
    G             % trained ft_mv_glmnet objects

    method = 'glmnet' % 'pls' (non-sparse), 'native' (sparse plus coupling), 'glmnet' (sparse based on glmnet package)
    
    alpha = 0.9;      % ridge versus lasso mixing parameter (alpha is a coupling matrix for native)
    lambda_min = 0.05; % regularization parameter (fraction of maximal lambda)
    nlambda = 20;     % number of steps to traverse the regularization path
    
  end

  methods
        
    function obj = ft_mv_sopls(varargin)
      
      obj = obj@ft_mv_predictor(varargin{:});
      
    end
    
    function obj = train(obj,X,Y)

      opts.alpha = obj.alpha;
      opts.lambda_min = obj.lambda_min;
      opts.nlambda = obj.nlambda;
      opts.verbose = obj.verbose;
      opts.method = obj.method;
        
      if ~strcmp(obj.method,'pls')
        
        [obj.A,obj.B,bias,obj.G] = ospls(X',Y',obj.nhidden,opts);
      
        % total bias;
        obj.C = bias * obj.A';
        
      else
        
        [XL,YL,XS,YS,beta,pctvar,mse,stats] = plsregress(X,Y,obj.nhidden);
        
        obj.B = stats.W;
        obj.A = YL;        
        obj.C = mean(Y,1) - mean(X,1)*obj.B*obj.A';
            
      end
      
    end
    
    function Y = test(obj,X)

      Y = bsxfun(@plus,X * obj.B(:,1:obj.nhidden) * obj.A(:,1:obj.nhidden)',obj.C);
       
    end
    
  end
  
end
