classdef ft_mv_glmnet < ft_mv_predictor
% FT_MV_GLMNET wrapper class for Friedman's glmnet package
%
% data is standardized inside the algorithm
%
% If obj.lambda is empty then a whole regularization path is estimated.
% Optimal performance is computed using a crossvalidator.
%
% 'family' can be 'gaussian' (linear regression), 'binomial' or 'multinomial' (logistic regression)
%
% Please use help glmnetSet to determine the possible options. Some
% examples:
%
% alpha   = 1;        % mixing parameter; alpha = 1 => L1, alpha < 1 => elastic net
% lambda  = [];       % L1 parameter; empty returns a whole path depending on nlambda and lambda_min
% nlambda = 100;      % number of models to evaluate in the regularization path
% lambda_min = 0.05;  % smallest value of lambda as a fraction of lambda_max
% standardize = 1;    % can be set to false
%
% EXAMPLE:
%
% [a,b,c,d] = ft_mv_test('mva',{ft_mv_glmnet('validator',ft_mv_crossvalidator('nfolds',5,'metric','logprob'))})
%
%   Copyright (c) 2010, Marcel van Gerven

  properties
    
    family = 'binomial' % 'gaussian' (linear regression), 'binomial' or 'multinomial' (logistic regression)
    
    validator           % crossvalidator object if a whole path is specified    
    performance         % performance results for the crossvalidator
   
    % replicated glmnetSet options (see help glmnetSet)
    weights
    alpha = 0.99
    nlambda = 100
    lambda_min = 0
    lambda = [];
    standardize = 1
    thresh = 1.0000e-04
    dfmax = 0
    pmax = 0
    exclude = []
    penalty_factor = []
    maxit = 100
    HessianExact = 0
    type = 'covariance'
    
    deregularize = 0; % retrain using unregularized regresson (not advised)
 
  end

  methods
    
    function obj = ft_mv_glmnet(varargin)
      
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
     
      % convert logical data to double (otherwise we get numerical problems)
      if islogical(X), X = double(X); end
      
      % multiple outputs
      if size(Y,2) > 1
        obj = ft_mv_noutput('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      if strcmp(obj.family,'gaussian')
        nclasses = 1;
      else
        nclasses = max(2,max(Y(:)));
        if nclasses > 2
          obj.family = 'multinomial';
        else
          obj.family = 'binomial';
        end
      end
      
      opts = [];
      opts.weights = obj.weights;
      opts.alpha = obj.alpha;
      opts.nlambda = obj.nlambda;
      opts.lambda_min = obj.lambda_min;
      opts.lambda = obj.lambda;
      opts.standardize = obj.standardize;
      opts.thresh = obj.thresh;
      opts.dfmax = obj.dfmax;
      opts.pmax = obj.pmax;
      opts.exclude = obj.exclude;
      opts.penalty_factor = obj.penalty_factor;
      opts.maxit = obj.maxit;
      opts.HessianExact = obj.HessianExact;
      opts.type = obj.type;
     
      if isempty(opts.lambda) && ~isempty(obj.validator)

        % use dummy classifier to determine the lambda path for all folds
        % otherwise we get different lambdas for different folds which
        % makes everything hard to compare
        dum = obj;
        dum.validator = [];
        dum = dum.train(X,Y);
        obj.lambda = dum.lambda;
        
        obj.validator.mva = obj;
        obj.validator.mva.validator = [];
        
        obj.validator = obj.validator.train(X,Y);
        
        % sometimes lambda is cut off by glmnet; this ensures the solutions
        % exist
        nsolutions = min(cellfun(@(x)(size(x,2)),obj.validator.post)) / nclasses;
        
        obj.performance = zeros(length(obj.validator.design),nsolutions);

        for j=1:nsolutions
          
          post = cellfun(@(x)(x(:,(j-1)*nclasses + (1:nclasses))),obj.validator.post,'UniformOutput',false);
          
          perf = ft_mv_performance(obj.validator.design,post,obj.validator.metric);
          
          if iscell(perf), perf = cell2mat(perf); end
          
          obj.performance(:,j) = perf;
          
        end
        
       % find best lambda
       [a,b] = max(obj.performance,[],2);

       % create lambda path with the best lambda (mean(lbest)) at the end
       % in order to ensure proper convergence to the correct solution
       lbest = zeros(1,length(b));
       lmax = zeros(1,length(b));
       for j=1:length(b)
         lbest(j) = obj.validator.mva{j}.mvmethods{end}.lambda(b(j));
         lmax(j) = obj.validator.mva{j}.mvmethods{end}.lambda(1);
       end
       opts.lambda = linspace(max(lmax),mean(lbest),50);
       obj.lambda = opts.lambda;
       
       u = unique(Y);
       
       if length(u) == 1
          % degenerate case
          
          if u==1
            obj.weights = [zeros(size(X,2),1); -1];
          else % u==2
            obj.weights = [zeros(size(X,2),1); 1];
          end
          
          if isempty(opts.lambda)
            obj.weights = repmat(obj.weights,[1 opts.nlambda]);
          end
          
          res.lambda = nan;
          
       else
         
         res = glmnet(X,Y,obj.family,opts);
         obj.weights = [res.beta(:,end); res.a0(end)]; 
       
       end
       
      else
      
        u = unique(Y);
        
        if length(u) == 1
          % degenerate case
          
          if u==1
            obj.weights = [zeros(size(X,2),1); -0.1];
          else % u==2
            obj.weights = [zeros(size(X,2),1); 0.1];
          end
          
          if isempty(opts.lambda)
            obj.weights = repmat(obj.weights,[1 opts.nlambda]);
          end
           
          res.lambda = nan;
        else          
          res = glmnet(X,Y,obj.family,opts);
          obj.weights = [res.beta; res.a0(:)'];
        end
        
      end
      
      obj.lambda  = res.lambda;
        
      if obj.deregularize
        W = obj.weights;
        for k=1:size(W,2)
          widx = W(1:(end-1),k) ~= 0;
          if any(widx)
            Xt = X(:,widx);
            if strcmp(obj.family,'gaussian')
              W([widx; true],k) = regress(Y,[Xt ones(size(Xt,1),1)]); % X \ Y
            else
              W([widx; true],k) = -logist2(Y-1,[Xt ones(size(Xt,1),1)]);
            end
          end
        end
        obj.weights = W;
      end
      
    end
    
    function Y = test(obj,X)

      X = [X ones(size(X,1),1)];
      
      if strcmp('gaussian',obj.family)

        Y = X * obj.weights;
        
      else
        
        if size(obj.weights,2) > 1 % regularization path
          
          Y = [];
          for j=1:size(obj.weights,2)
            Z = exp(X * [obj.weights(:,j) zeros(size(X,2),1)]);
            Z =  1 - bsxfun(@rdivide,Z,sum(Z,2));
            Y = [Y Z];
          end
          
        else
          Y = exp(X * [obj.weights zeros(size(X,2),1)]);
          Y = 1 - bsxfun(@rdivide,Y,sum(Y,2));
        end
        
        % MULTINOMIAL NOT (YET) SUPPORTED
      
      end
      
    end
    
    function [m,desc] = model(obj)
      % return the parameters wrt a class label in some shape
      
      % weight-vector for class 1 is always zero
      m =  mat2cell(obj.weights',ones(1,size(obj.weights,2)),size(obj.weights,1));
      m{length(m)+1,1} = zeros(size(m{1}));
      
      desc = cell(length(m),1);
      for j=1:length(m)

        m{j} = m{j}(:,1:(size(m{j},2)-1));
        
        desc{j} = sprintf('regression coefficients for class %d',j);
      end
      
    end
    
  end
end
