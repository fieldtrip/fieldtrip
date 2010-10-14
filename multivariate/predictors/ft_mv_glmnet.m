classdef ft_mv_glmnet < ft_mv_predictor
% FT_MV_GLMNET wrapper class for Friedman's glmnet package
%
% data is standardized inside the algorithm
%
% If obj.lambda is empty then a whole regularization path is estimated.
% Optimal performance is computed using a crossvalidator.
%
% EXAMPLE:
%
% [a,b,c,d] = ft_mv_test('mva',{ft_mv_glmnet('cv',ft_mv_crossvalidator('nfolds',5,'metric','accuracy'))})
%
%
%   Copyright (c) 2010, Marcel van Gerven

  properties
    
    alpha   = 1;        % mixing parameter; alpha = 1 => L1, alpha < 1 => elastic net
    
    lambda  = [];       % L1 parameter; empty returns a whole path depending on nlambda and lambda_min
  
    nlambda = 100;      % number of models to evaluate in the regularization path
    lambda_min = 0.05;  % smallest value of lambda as a fraction of lambda_max
    
    type = 'logistic'   % 'linear' or 'logistic' regression
    
    weights             % regression coefficients
    
    cv                  % crossvalidator object if a whole path is specified    
    performance         % performance results for the crossvalidator
    
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
     
      % multiple outputs
      if size(Y,2) > 1
        obj = ft_mv_noutput('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      if strcmp(obj.type,'linear')
        family = 'gaussian';
        nclasses = 1;
      else
        nclasses = max(2,max(Y(:)));
        if nclasses > 2
          family = 'multinomial';
        else
          family = 'binomial';
        end
      end
      
      opts = glmnetSet;
      
      opts.alpha = obj.alpha; % mixing parameter
      
      if ~isempty(obj.lambda), opts.lambda = obj.lambda; end 

      if isempty(obj.lambda) && ~isempty(obj.cv)

        obj.cv.mva = obj;
        obj.cv.mva.cv = [];
        
        obj.cv = obj.cv.train(X,Y);
        
        nsolutions = size(obj.cv.post{1},2);
        
        obj.performance = zeros(length(obj.cv.design),nsolutions/nclasses);

        for j=1:(nsolutions/nclasses)
          
          post = cellfun(@(x)(x(:,(j-1)*nclasses + (1:nclasses))),obj.cv.post,'UniformOutput',false);
          
          perf = ft_mv_performance(obj.cv.design,post,obj.cv.metric);
          
          if iscell(perf), perf = cell2mat(perf); end
          
          obj.performance(:,j) = perf;
          
        end
        
       % find best lambda
       [a,b] = max(obj.performance,[],2);
       
       lbest = zeros(1,length(b));
       lmax = zeros(1,length(b));
       for j=1:length(b)
         lbest(j) = obj.cv.mva{j}.mvmethods{end}.lambda(b(j));
         lmax(j) = obj.cv.mva{j}.mvmethods{end}.lambda(1);
       end
       
       
       opts.lambda = linspace(mean(lmax),mean(lbest),20);
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
         
         res = glmnet(X,Y,family,opts);
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
          res = glmnet(X,Y,family,opts);
          obj.weights = [res.beta; res.a0(:)'];
        end
        
      end
      
      obj.lambda  = res.lambda;
          
    end
    
    function Y = test(obj,X)

      X = [X ones(size(X,1),1)];
      
      if strcmp('linear',obj.type)

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
        
        % MULTINOMIAL NOT SUPPORTED
      
      end
      
    end
    
    function [m,desc] = model(obj)
      % return the parameters wrt a class label in some shape
      
      % weight-vector for class 1 is always zero
      m =  mat2cell(obj.weights',ones(1,size(obj.weights,2)),size(obj.weights,1));
      m{length(m)+1,1} = zeros(size(m{1}));
      
      desc = cell(length(m),1);
      for j=1:length(m)

        % RETURNS ABSOLUTE VALUES
        m{j} = abs(m{j}(:,1:(size(m{j},2)-1)));
        
        desc{j} = sprintf('magnitude of the regression coefficients for class %d',j);
      end
      
    end
    
  end
end
