classdef graphnet < dml.method
% GRAPHNET native implementation of graphnet algorithm.
%
%   DESCRIPTION
%   Elastic net linear and logistic regression. Note that this algorithm
%   allows for the GraphNet generalization that allows coupling between
%   features. Requires specification of the lasso penalty term L1 and
%   (optionally) the ridge penalty term L2. The latter can also be a matrix
%   to allow coupling of variables.
%
%   REFERENCE Regularization paths for generalized linear models via
%   coordinate descent by Friedman et al.
%
%   Grosenick L, Klingenberg B, Knutson B. A family of interpretable
%   multivariate models for regression and classification of whole-brain
%   fMRI data. stanford.edu.
%
%   EXAMPLE 
%   X = rand(10,20); Y = [1 1 1 1 1 2 2 2 2 2]'; 
%   m = dml.graphnet('family','binomial','L1',0.1);
%   m = m.train(X,Y); 
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
   weights % regression weights (offset last)
   
   family = 'gaussian' % gaussian, binomial, or multinomial
   
   L1 % lasso penalty
   L2 = 1e-6 % nfeatures x nfeatures ridge penalty    
      
   maxiter = 1e4; % maximum number of iterations for native elastic net
   tolerance = 1e-4; % tolerance in the error for native elastic net
      
   conv % plot of the convergence of the parameters sum(abs(beta-betaold));
   
   df % degrees of freedom
   
  end
  
  methods
    
    function obj = graphnet(varargin)

      obj = obj@dml.method(varargin{:});

    end
    
    function obj = train(obj,X,Y)
        
      % handle multiple datasets
      if iscell(X)
        obj = dml.ndata('method',obj);
        obj = obj.train(X,Y);
        return;
      end
        
      if isempty(obj.L1), error('unspecified L1 parameter'); end
      
      % transform scalar L2 to matrix
      if isscalar(obj.L2), obj.L2 = obj.L2*eye(size(X,2)); end
      
      % options for elastic net
      opt = [];
      opt.offset = 1;
      opt.maxiter = obj.maxiter;
      opt.tol = obj.tolerance;
      
      switch obj.family
        
        case 'gaussian'
          
          if obj.L1 == 0 % unregularized or ridge regression
            
            l2 = obj.L2; l2(end+1,end+1) = 0;
            X = [X ones(size(X,1),1)];
            R = chol(X'*X + l2);
            obj.weights = R\(R'\(X'*Y));
            
          else % elastic net linear regression
            
            if obj.restart || isempty(obj.weights)
              [beta,beta0,c] = elasticlin(X',Y',obj.L1,obj.L2,opt);
              obj.conv = c;
            else
              % run model starting at the previous weights
              [beta,beta0,c] = elasticlin(X',Y',obj.L1,obj.L2,opt,obj.weights(1:(end-1))',obj.weights(end));
              obj.conv = c;
            end
            obj.weights = [beta; beta0]';
            
          end
          
          idx = (obj.weights ~= 0); idx(end)=0;
          obj.df = trace(X(:,idx)*inv(X(:,idx)'*X(:,idx) + obj.L2(idx,idx))*X(:,idx)');
          
        case 'binomial'
          
          if obj.restart || isempty(obj.weights)
            [beta,beta0,c] = elasticlog(X',Y' - min(Y),obj.L1,obj.L2,opt);
            obj.conv = c;
          else
            % run model starting at the previous weights
            [beta,beta0,c] = elasticlog(X',Y' - min(Y),obj.L1,obj.L2,opt,obj.weights(1:(end-1))',obj.weights(end));
            obj.conv = c;
          end
          obj.weights = [beta; beta0]';
                    
        otherwise
          
          error('unrecognized type');
          
      end
      
    end
    
    function Y = test(obj,X)

      switch obj.family
        
        case 'gaussian'
          
            Y = [X ones(size(X,1),1)] * obj.weights';
          
        case 'binomial'
          
          Y = (1./(1+ exp([X ones(size(X,1),1)] * obj.weights')));
          Y = [Y 1-Y];
          
        otherwise
          error('unrecognized kind');
      
      end
      
    end

    function m = model(obj)
    % returns 
    %
    % m.weights regression coefficients
    
      m.weights = obj.weights(1:(end-1))';
      m.bias = obj.weights(end);
      
    end
  
  end
  
  methods(Static=true)
    
    function p = lambdapath(X,Y,family,nsteps,lmin)
    
      if nargin < 4, nsteps = 50; end
      if nargin < 5, lmin = 1e-4; end
      
      switch family
        case 'gaussian'
          lmax = max(abs(X' * (Y - mean(Y)))) / size(X,1); % subtract offset
        case 'binomial'
          lmax = max(abs(X' * (Y-1.5)));
        otherwise
          error('unrecognized type');
      end
      p = exp(linspace(log(lmax),log(lmin*lmax),nsteps));
      
    end
    
    function K = laplacian(dims,s)
      % return matrix Laplacian for a multidimensional array of size dims and strength s
    
      nfeatures = prod(dims);
      
      cdim = cumprod(dims);
      tdim = [1 cdim(1:(end-1))];
      
      K = spalloc(nfeatures,nfeatures,nfeatures+nfeatures*ceil(2^length(dims)/2));
      
      for i=1:nfeatures
        
        for d=1:length(dims)
          
          % get neighbouring features in all dimensions
          
          % ignore the end boundary
          if mod(ceil(i/tdim(d)),dims(d)) ~= 0
            
            nbr = i + tdim(d);
            if nbr <= nfeatures
              K(i,nbr) = -1;
            end
            
          end
          
        end
      end
      
      K = K + K';
      
      K(1:(nfeatures+1):numel(K)) = -sum(K,2);

      K = s*K;
      
    end
      
  end
  
end
