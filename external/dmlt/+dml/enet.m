classdef enet < dml.method
% ENET efficient elastic net algorithm.
%
%   DESCRIPTION
%   Fast implementation of elastic net linear and logistic regression by
%   Friedman et al. Without arguments, this method automatically determines 
%   the lambda path and uses cross-validation to find the best lambda. The
%   parameter 'family' determines whether to use linear regression
%   (gaussian) or logistic regression (binomial or multinomial). The mixing
%   parameter alpha controls the influence of the 11 (alpha=1) and l2 regularizer
%   (alpha = 0). Negative weights indicate the class with label 1. Positive
%   weights indicate the class with label 2.
%
%   REFERENCE
%   Regularization paths for generalized linear models via coordinate descent 
%   by Friedman et al.
%
%   EXAMPLE
%   X = randn(100,20); Y = [ones(50,1); 2*ones(50,1)]; 
%   X(1:50,1:10) = X(1:50,1:10)+1; X(51:100,11:20) = X(51:100,11:20)+1;
%   m = dml.enet('family','binomial');
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   X = rand(15,20); Y = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]'; X(1:5,:) = X(1:5,:)+1;
%   m = dml.enet('family','multinomial');
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
   weights % regression weights (offset last)
   
   family = 'gaussian' % gaussian, binomial, or multinomial
   
   lambda;      % used lambda(s)
   alpha = 0.5;   % glmnet mixing parameter
   validator = dml.crossvalidator; % cross-validator

   performance % cross-validated decoding performance for elements on lambda path
   path % lambda path followed when estimating decoding performance
   
   df % degrees of freedom (not computed if zero)
   
  end
  
  methods
    
    function obj = enet(varargin)

      obj = obj@dml.method(varargin{:});

    end
    
    function obj = train(obj,X,Y)
      
      % handle multiple datasets
      if iscell(X)
        obj = dml.ndata('method',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % missing data
      if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
      
       % multiple outputs
      if size(Y,2) > 1
        obj = dml.noutput('method',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % convert data to double (otherwise we get numerical problems)
      if ~isa(X(1),'double'), X = double(X); end
      if ~isa(Y(1),'double'), Y = double(Y); end
      
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
          
        elseif obj.alpha == 0 && strcmp(obj.family,'gaussian')
          
          lambdas = [obj.lambda*ones(size(X,2),1); 0];
          X = [X ones(size(X,1),1)];
          R = chol(X'*X + diag(lambdas));
          obj.weights = R\(R'\(X'*Y));
          
        elseif obj.lambda == 0 && strcmp(obj.family,'binomial')
          
          obj.weights = logist2(Y-1,[X ones(size(X,1),1)]);
          
        elseif obj.lambda == 0 && strcmp(obj.family,'multinomial')
          
          P = zeros(length(Y),max(Y));
          for i=1:length(Y)
            P(i,Y(i)) = 1;
          end
          
          obj.weights = logistK(P',[X ones(size(X,1),1)]')';
          %obj.weights = mnrfit(X,P,'model','nominal','interactions','on');
          
        else
          
            opts = glmnetSet;
            if ~isempty(obj.alpha), opts.alpha = obj.alpha; end
            if ~isempty(obj.lambda), opts.lambda = obj.lambda; end
          
            obj.weights = obj.estimate(X,Y,obj.family,opts);
            if strcmp(obj.family,'multinomial')
              if ndims(obj.weights)==3 && size(obj.weights,3)>1
                obj.weights = obj.weights(:,:,end);
              end
            else
              if ndims(obj.weights)==2 && size(obj.weights,2)>1
                obj.weights = obj.weights(:,end);
              end
            end
            
        end
        
       else
        
        opts = glmnetSet;
        if ~isempty(obj.alpha), opts.alpha = obj.alpha; end
        if ~isempty(obj.lambda), opts.lambda = obj.lambda; end
        
        if ~isscalar(obj.lambda) && ~isempty(obj.validator)
          % find best lambda by using the cross-validator
          
          % predetermine the lambda path for all folds
          % otherwise we get different lambdas for different folds which
          % makes everything hard to compare
          opts.iter = 0;
          res = glmnet(X,Y,obj.family,opts);
          obj.lambda = res.lambda;
          
          if isempty(obj.validator.stat)
            if strcmp(obj.family,'gaussian')
              obj.validator.stat='-RMS';
            else
              obj.validator.stat='accuracy';
            end
          end
          
          obj.validator.mva = obj;
          obj.validator.mva.validator = [];
          obj.validator = obj.validator.train(X,Y);
          
          % sometimes lambda is cut off by glmnet; this ensures the solutions exist
          nsolutions = min(cellfun(@(x)(size(x,2)),obj.validator.result)) / nclasses;
          
          obj.performance = zeros(length(obj.validator.design),nsolutions);
          
          for j=1:nsolutions
            
            post = cellfun(@(x)(x(:,(j-1)*nclasses + (1:nclasses))),obj.validator.result,'UniformOutput',false);
            
            perf = dml.statistic(obj.validator.stat,cell2mat(obj.validator.design),cell2mat(post));
            
            if iscell(perf), perf = cell2mat(perf); end
            
            obj.performance(:,j) = perf;
            
          end
          obj.path  = obj.lambda(1:nsolutions);
          
          % find best lambda; smallest lambda whose performance is within one standard
          % error from the best performance (Hastie's one standard error rule)
          %[a,b] = max(obj.performance,[],2);
          if size(obj.performance,1)==1
            [a,b] = max(obj.performance);
          else
            mp = mean(obj.performance);
            [maxp,midx] = max(mp);
            ep = std(obj.performance) ./ sqrt(size(obj.performance,1));
            b = find(mp >= (maxp - ep(midx)),1,'first');
            if isempty(b), b=1; end
          end
          
          % create lambda path with the best lambda (mean(lbest)) at the end
          % in order to ensure proper convergence to the correct solution
          lbest = zeros(1,length(b));
          lmax = zeros(1,length(b));
          for j=1:length(b)
            lbest(j) = obj.validator.mva{j}.lambda(b(j));
            lmax(j) = obj.validator.mva{j}.lambda(1);
          end
          obj.lambda = mean(lbest);
          opts.lambda = linspace(max(lmax),mean(lbest),50);
          
          % create unique path (1 element if lbest=lmax)
          opts.lambda = sort(unique(opts.lambda),'descend');
          
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
            
          else
            
            obj.weights = obj.estimate(X,Y,obj.family,opts);
            
          end
          
          if strcmp(obj.family,'multinomial')
            obj.weights = obj.weights(:,:,end);
          else
            obj.weights = obj.weights(:,end);
          end
          obj.performance = mean(obj.performance,1);
          
        else % return result for specified lambda; no cross-validation
          
          if isscalar(obj.lambda)
            
            obj.weights = obj.estimate(X,Y,obj.family,opts);
            if strcmp(obj.family,'multinomial')
              if ndims(obj.weights)==3 && size(obj.weights,3)>1
                obj.weights = obj.weights(:,:,end);
              end
            else
              if ndims(obj.weights)==2 && size(obj.weights,2)>1
                obj.weights = obj.weights(:,end);
              end
            end
            
          else
            
            [obj.weights,obj.path] = obj.estimate(X,Y,obj.family,opts);
            obj.lambda = obj.path;
            
          end         
          
        end
        
      end
      
      % compute degrees of freedom for Gaussian case
      if isempty(obj.df) && size(obj.weights,2)==1 && strcmp(obj.family,'gaussian')
        
        idx = (obj.weights ~= 0); idx(end)=0;
        obj.df = trace(X(:,idx)*inv(X(:,idx)'*X(:,idx) + obj.lambda*(1-obj.alpha)*eye(sum(idx)))*X(:,idx)');

      end

               
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
    
    function m = model(obj)
      % returns
      %
      % m.weights regression coefficients
      % m.bias bias term
      % m.lambda selected lambda value
      
      m.weights = obj.weights(1:(end-1),:);
      m.bias = obj.weights(end,:);
      m.lambda = obj.lambda;
      
    end
    
  end
  
  methods(Access=private)
    
    function [weights,lambda] = estimate(obj,X,Y,family,opts)
      % safeguarded estimation of glmnet
      
      try
        
        res = glmnet(X,Y,obj.family,opts);
        
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
      
      if strcmp(family,'multinomial')
        x = cell2mat(cellfun(@(x)(reshape(x,[1 size(res.beta{1})])),res.beta,'UniformOutput',false)');
        weights = permute(cat(2,x,reshape(res.a0,[size(x,1) 1 size(x,3)])),[2 1 3]);
      else
        weights = [res.beta; res.a0(:)'];
      end
      
      lambda  = res.lambda;

    end
    
  end    
  
end
