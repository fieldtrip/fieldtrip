classdef ft_mv_glmnet < ft_mv_predictor
% FT_MV_GLMNET wrapper class for Friedman's glmnet package
%
% data is standardized inside the algorithm and bias term is handled
% automatically.
%
% If obj.lambda is empty then a whole regularization path is estimated. 
% Optimal performance can be computed using a crossvalidator.
%
% If obj.lambda is non-empty then the solution for obj.lambda(end) is returned.
% Specifying a decreasing path is more stable than specifying one solution.
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
% [a,b,c] = ft_mv_test('mva',{ft_mv_glmnet('validator',ft_mv_crossvalidator('nfolds',5,'metric','logprob'))})
%
% NOTE:
% In case you encounter a failure because libgfortran is not installed,
% please follow the instructions from here <http://www.cs.ubc.ca/~hoffmanm/matlab.html>
% to download the libgfortran 4.2 tarball from here
% <http://r.research.att.com/tools/> listed under OS X 10.6 and copy libgfortran* 
% files to /usr/local/lib.
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
    lambda_min = 0  % minimal lambda as fraction of lambda_max
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
    
    % retrain using unregularized regression (not advised; should be
    % handled by weighted regression)
    deregularize = 0; 
 
    adaptive = false;
    
    % 'glmnet' or 'native' (latter only for 'gaussian' and with limited
    % parameters. 'native' is useful when we want lambda and alpha to
    % represent the L1 and L2 penalty separately and if we want alpha to
    % implement spatial smoothing. alpha = 0 then means no L2 regularization
    % opposite to the normal case. Also requires smaller lambda values in
    % the lambda path. Currently, using 'glmnet' is recommended.
    method = 'glmnet'; 
    
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
     
      % convert data to double (otherwise we get numerical problems)
      if ~isa(X(1),'double'), X = double(X); end
      if ~isa(Y(1),'double'), Y = double(Y); end
      
      % multiple outputs
      if size(Y,2) > 1
        obj = ft_mv_noutput('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      if strcmp(obj.family,'gaussian')
        nclasses = 1;
      else
        if any(rem(Y(:),1))
          error('expecting class labels for logistic regression');
        end
        nclasses = max(2,max(Y(:)));
        if nclasses > 2
          obj.family = 'multinomial';
        else
          obj.family = 'binomial';
        end
      end
      
      % handle some special cases
      if isscalar(obj.lambda)
        
        if obj.lambda == 0 && strcmp(obj.family,'gaussian')
        
          obj.weights = regress(Y,[X ones(size(X,1),1)]); % X \ Y
          
          return
        
        elseif obj.alpha == 0 && strcmp(obj.family,'gaussian')
          
          lambdas = [obj.lambda*ones(size(X,2),1); 0];
          X = [X ones(size(X,1),1)];
          R = chol(X'*X + diag(lambdas));
          obj.weights = R\(R'\(X'*Y));
          
          return
          
        elseif obj.lambda == 0 && strcmp(obj.family,'binomial')
          
          obj.weights = -logist2(Y-1,[X ones(size(X,1),1)]);
          
          return
          
        end
        
      end
      
      if obj.adaptive 
        % adaptive lasso; implemented by using univariate ols regression
        % estimates as penalty weights
        
        if obj.verbose
          fprintf('estimating penalty factor for adaptive method\n');
        end
        
        obj.penalty_factor = zeros(1,size(X,2));
        
         if strcmp(obj.family,'gaussian')
           
           for j=1:size(X,2)
             w = regress(Y,[X(:,j) ones(size(X,1),1)]); % X \ Y
             obj.penalty_factor(j) = 1/abs(w(1))^obj.adaptive;
           end     
           
         elseif strcmp(obj.family,'binomial')
           
           for j=1:size(X,2)
             w = -logist2(Y-1,[X(:,j) ones(size(X,1),1)]);
             obj.penalty_factor(j) = 1/abs(w(1))^obj.adaptive;
           end
           
         else % multinomial

           YY = zeros(size(X,1),nclasses);
           for j=1:size(X,1)
             YY(j,Y(j)) = 1;
           end
           for j=1:size(X,2)
             w = -logistK(YY',[X(:,j) ones(size(X,1),1)]');
             obj.penalty_factor(j) = 1/mean(abs(w(:,1)))^obj.adaptive; % mean absolute coefficient
           end
           
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
      
      if strcmp(obj.method,'native') && isscalar(opts.alpha)
        opts.alpha = opts.alpha*diag(ones(size(X,2),1));
      end

      if isempty(opts.lambda) && ~isempty(obj.validator)

        % use dummy classifier to determine the lambda path for all folds
        % otherwise we get different lambdas for different folds which
        % makes everything hard to compare
        
        if obj.verbose, fprintf('determining lambda path\n'); end
        
        % explicit computation of lambda_path by glmnet
%         dum = obj;
%         dum.validator = [];
%         dum = dum.train(X,Y);
%         obj.lambda = dum.lambda;

        obj.lambda = obj.lambda_path(X,Y)';
        
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
        
       % find best lambda; smallest lambda whose performance is within one standard
       % error from the best performance (Hastie's one standard error rule)
       %[a,b] = max(obj.performance,[],2);
       if size(obj.performance,1)==1
         [a,b] = max(obj.performance);
       else
         mp = mean(obj.performance);
         [maxp,midx] = max(mp);
         ep = std(obj.performance) ./ size(obj.performance,1);
         b = find(mp >= (maxp - ep(midx)),1,'first');
         if isempty(b), b=1; end
       end
       
       % create lambda path with the best lambda (mean(lbest)) at the end
       % in order to ensure proper convergence to the correct solution
       lbest = zeros(1,length(b));
       lmax = zeros(1,length(b));
       for j=1:length(b)
         lbest(j) = obj.validator.mva{j}.mvmethods{end}.lambda(b(j));
         lmax(j) = obj.validator.mva{j}.mvmethods{end}.lambda(1);
       end
       opts.lambda = linspace(max(lmax),mean(lbest),50);
       
       % create unique path (1 element if lbest=lmax)
       opts.lambda = sort(unique(opts.lambda),'descend');
       obj.lambda  = opts.lambda;
       
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
         
         try
           
           if strcmp(obj.method,'glmnet')
           
             res = glmnet(X,Y,obj.family,opts);

           else % native implementation
             
             % implemented for gaussian only at the moment
             assert(strcmp(obj.family,'gaussian')); 
             
             beta = zeros(size(X,2),1);
             beta0 = 0;
             for j=1:numel(opts.lambda)
               if obj.verbose
                 fprintf('estimating model %d of %d with lambda = %.2f\n',j,numel(opts.lambda),opts.lambda(j));
               end
               [beta,beta0] = elastic(X',Y',opts.lambda(j),opts.alpha,opts,beta,beta0);
               obj.weights = cat(2,obj.weights,[beta; beta0]);
             end
             res.lambda = opts.lambda;
             res.beta = obj.weights(1:(end-1),:);
             res.a0 = obj.weights(end,:);

             
           end
           
         
         catch
           
           % probably returned the empty model
           warning(lasterr);
           
           res.lambda = obj.lambda(1);
           nclasses = max(Y);
           if nclasses>2
             res.beta = cell(1,nclasses);
             for c=1:nclasses
               res.beta{c} = zeros(size(X,2),1);
               res.a0 = zeros(nclasses,1);
             end
           else
             res.beta = zeros(size(X,2),1);
             res.a0   = 0;
           end
           
         end
         
         if strcmp(obj.family,'multinomial')
            obj.weights = zeros(size(res.beta{1},1)+1,size(res.a0,1));
            for c=1:size(res.a0,1)
              obj.weights(:,c) = [res.beta{c}(:,end); res.a0(c,end)];
            end
         else
            obj.weights = [res.beta(:,end); res.a0(end)]; 
          end
       
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
          
          try
          
            if strcmp(obj.method,'glmnet')

              res = glmnet(X,Y,obj.family,opts);

            else % native implementation
              
              % implemented for gaussian only at the moment
              assert(strcmp(obj.family,'gaussian'));
              
              beta = zeros(size(X,2),1);
              beta0 = 0;
              for j=1:numel(opts.lambda)
                if obj.verbose
                 fprintf('estimating model %d of %d with lambda = %.2f\n',j,numel(opts.lambda),opts.lambda(j));
                end
                [beta,beta0] = elastic(X',Y',opts.lambda(j),opts.alpha,opts,beta,beta0);
                obj.weights = cat(2,obj.weights,[beta; beta0]);
              end
              res.lambda = opts.lambda;
              res.beta = obj.weights(1:(end-1),:);
              res.a0 = obj.weights(end,:);
              
            end
              
            %  can't be used since the validator part requires the full
            %  path!
%             if strcmp(obj.family,'multinomial')
%               for c=1:length(res.beta)
%                 res.beta{c} = res.beta{c}(:,end);
%               end
%               res.a0 = res.a0(:,end);
%             else
%               res.beta = res.beta(:,end);
%               res.a0 = res.a0(end);
%             end
            
          catch
            
            warning(lasterr);
            
            % probably returned the empty model
            res.lambda = obj.lambda(1);
            nclasses = max(Y);
            if strcmp(obj.family,'multinomial')
              res.beta = cell(1,nclasses);
              for c=1:nclasses
                res.beta{c} = zeros(size(X,2),1);
                res.a0 = zeros(nclasses,1);
              end
            else
              res.beta = zeros(size(X,2),1);
              res.a0   = 0;
            end
            
          end
          
          if strcmp(obj.family,'multinomial')
            x = cell2mat(cellfun(@(x)(reshape(x,[1 size(res.beta{1})])),res.beta,'UniformOutput',false)');
            obj.weights = permute(cat(2,x,reshape(res.a0,[size(x,1) 1 size(x,3)])),[2 1 3]);
          else
            obj.weights = [res.beta; res.a0(:)'];
          end
        end
        
      end
      
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
      
      obj.lambda  = res.lambda;
              
    end
    
    function Y = test(obj,X)

      X = [X ones(size(X,1),1)];
      
      if strcmp(obj.family,'gaussian')

        Y = X * obj.weights;
        
      elseif strcmp(obj.family,'binomial')
        
        if size(obj.weights,2) > 1% regularization path
          
          Y = [];
          for j=1:size(obj.weights,ndims(obj.weights))
            Z = exp(X * [obj.weights(:,j) zeros(size(X,2),1)]);
            Z =  1 - bsxfun(@rdivide,Z,sum(Z,2));
            Y = [Y Z];
          end
          
        else
          Y = exp(X * [obj.weights zeros(size(X,2),1)]);
          Y = 1 - bsxfun(@rdivide,Y,sum(Y,2));
        end
        
      else
        
        if size(obj.weights,3) > 1% regularization path
          
          Y = [];
          for j=1:size(obj.weights,ndims(obj.weights))
            Z = exp(X * obj.weights(:,:,j));
            Z =  bsxfun(@rdivide,Z,sum(Z,2));
            Y = [Y Z];
          end
          
        else
          Y = exp(X * obj.weights);
          Y = bsxfun(@rdivide,Y,sum(Y,2));
        end
        
      end
      
    end
    
    function [m,desc] = model(obj)
      % return the parameters wrt a class label in some shape
      
      m = {};
      if ~strcmp(obj.family,'multinomial')
        if size(obj.weights,2)==1 % only return for one model and not for path
          m =  mat2cell(obj.weights',ones(1,size(obj.weights,2)),size(obj.weights,1));
          % weight-vector for class 1 is always zero
          m{length(m)+1,1} = zeros(size(m{1}));
        end
      else
        nclasses = size(obj.weights,2);
        m = cell(nclasses,1);
        for c=1:nclasses
          m{c} = obj.weights(:,c)';
        end
      end
      
      desc = cell(length(m),1);
      for j=1:length(m)

        m{j} = m{j}(:,1:(size(m{j},2)-1));
        
        desc{j} = sprintf('regression coefficients for class %d',j);
      end
      
    end
    
    function p = lambda_path(obj,X,Y)
      % compute lambda values that are used when computing the lambda path
      % only valid for binomial and gaussian case
      
      sz = size(X);
      
      assert(numel(sz)==2);
      
      N = sz(1);
      M = prod(sz(2:end));
      
      epsilon = obj.lambda_min;
      if isempty(epsilon) || epsilon==0
        if N > M
          epsilon = 0.0001;
        else
          epsilon = 0.05;
        end
      end
      
      X = zscore(X);
      Y = zscore(Y);
      
      if strcmp(obj.method,'glmnet')
        lmax = max(abs(sum(bsxfun(@times,X,Y))))./(obj.alpha*N);
      else
        lmax = max(abs(sum(bsxfun(@times,X,Y))))/N; % choose solution independent of alpha
      end
      
      p = exp(linspace(log(lmax),log(epsilon*lmax),obj.nlambda));
      
    end
    
    function b = bias(obj)
      % return bias term; required function for ft_mv_sopls
      b = obj.weights(end);
    end
    
  end
  
end
