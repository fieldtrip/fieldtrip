classdef blogreg < dml.method
% BLOGREG Bayesian logistic regression.
%
%   DESCRIPTION
%   Multivariate Laplace prior supports spatiotemporal interactions,
%   multitask and mixed effects models. The scale property specifies the
%   regularization. The bigger the scale, the less the regression
%   coefficients will be regularized towards zero. If multiple scales are
%   used then the optimal one will be selected based on the log model
%   evidence (recorded by the logp property). Note: a bias term is added
%   to the model and should not be included explicitly
%
%   EXAMPLE
%   In the following examples we work with the following input data:
%
%   rand('seed',1); randn('seed',1);
%   X1 = rand(10,5,10); X2 = X1 + 0.1*randn(size(X1));
%   Y1 = [1 1 1 1 1 2 2 2 2 2]'; Y2 = [1 1 1 1 2 1 2 2 2 2]';
%
%   Sometimes we assume the input data is just a region of interest that is
%   specified by some mask. E.g.
%
%   rand('seed',1); randn('seed',1);
%   mask = rand(5,10)>0.3;
%   X1m = rand(10,sum(mask(:))); X2m = X1m + 0.1*randn(size(X1m)); % create a subset of the data
%
%   creates data X1 (and X2) whose columns stand for those elements in the
%   mask that are equal to one. The representation for the output stays the
%   same.
%
%   In the examples use spy(f.prior) to look at the structure of the coupling
%   matrix
%
%   f = dml.blogreg('scale',logspace(-3,0,4));
%   f = f.train(X1,Y1);
%   f.test(X1)
%
%   Blogreg allows features to become coupled. This is done through the
%   coupling property, which specifies for each input dimension the coupling
%   strength in that dimension. E.g. coupling = [100 0 100] will strongly couple
%   neighboring features in the first and second input dimension. Also,
%   the input dimensions of the original data must be specified through the indims property.
%
%   f = dml.blogreg('indims',[5 10],'coupling',[100 100]);
%   f = f.train(X1,Y1);
%   f.test(X1)
%
%   In case the input data X is just a region of interest then we can use a
%   mask to indicate which features from the original volume of data are
%   represented by X. This still allows us to use the above approach to
%   specify the coupling.
%
%   f = dml.blogreg('indims',[5 10],'coupling',[100 100],'mask',mask);
%   f = f.train(X1m,Y1);
%   f.test(X1m)
%
%   In the following we discuss some more exotic uses of blogreg, which
%   is not required by the typical user.
%
%   Multitask learning (multitask = true) is implemented by augmenting the data matrix as
%   [ T1     0   ]
%   [ 0      T2  ]
%   etc and coupling the tasks through the taskcoupling property. Data X and Y must be given by a cell-array.
%
%   f = dml.blogreg('multitask',1);
%   f = f.train({X1 X2},{Y1 Y2});
%   f.test({X1 X2})
%
%   This can also be combined with coupling of the features themselves.
%
%   f = dml.blogreg('multitask',1,'indims',[5 10],'coupling',[100 100],'mask',mask);
%   f = f.train({X1m X2m},{Y1 Y2});
%   f.test({X1m X2m})
%
%   Blogreg also supports a mixed effects model (mixed = true)
%   of the form
%   [ M1 M1 0 ]
%   [ M2 0 M2 ]
%   where the first column contains the "fixed effects" and the remaining
%   columns the "random effects". The basic idea is that if prediction is
%   supported by a fixed effect then it will be chosen since this incurs the
%   smallest penalty in terms of sparseness. Data should be given by a
%   cell-array.
%
%   f = dml.blogreg('mixed',1);
%   f = f.train({X1 X2}',{Y1 Y2}');
%   f.test({X1 X2})
%
%   This can also be combined with coupling of the features themselves.
%   If a coupling is specified then this coupling will only operate on the
%   fixed effects using mixed=1 and on both the fixed and random effects
%   using mixed=2.
%
%   f = dml.blogreg('mixed',1,'indims',[5 10],'coupling',[100 100],'mask',mask);
%   f = f.train({X1m X2m}',{Y1 Y2}');
%   f.test({X1m X2m})
%
%   If the input data is a cell-array of cell-arrays then we assume a mixed
%   effects model for multiple tasks. The output will have the same
%   structure. E.g. model{i}{j,k} will be the j-th model for the k-th mixed
%   effect in the i-th subject. Note that k=1 is the fixed effect and k>1 are
%   the random effects.
%
%   f = dml.blogreg('mixed',1,'multitask',1);
%   f = f.train({{X1 X2} {X1 X2}},{{Y1 Y2} {Y1 Y2}});
%   f.test({{X1 X2} {X1 X2}})
%
%   If mixed = 1 then the only coupling  (spatial or multitask) will be for
%   the fixed effects part. For mixed = 2, also the random effects will be coupled.
%
%   f = dml.blogreg('mixed',1,'multitask',1,'coupling',[100 100],'indims',[5 10],'mask',mask);
%   f = f.train({ {X1m X2m} {X1m X2m} },{ {Y1 Y2} {Y1 Y2} });
%   f.test({{X1m X2m} {X1m X2m}})
%
%
%   REFERENCE
%
%   When using this method please refer to the following:
%
%   van Gerven et al. Efficient Bayesian multivariate fMRI analysis using a
%   sparsifying spatio-temporal prior. Neuroimage, 2010
%
%   van Gerven & Simanova. Concept classification with Bayesian multitask
%   learning, NAACL, 2010
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  
  properties
    
    % number of features; constant over tasks
    nfeatures
    
    % prior precision matrix of the auxiliary variables
    prior
    
    % scale of the bias term (bias term will be appended to the model)
    % note: internally scale and precbias (scale of bias term) will be
    % inverted. This is handled automatically but is a concern when
    % maually specifying and scaling the prior
    precbias = [];
    
    % (either 'probit' or gaussian 'quadrature' to approximate posterior)
    approximation = 'quadrature'
    
    % number of weights for gaussian quadrature
    nweights = 100;
    
    % coupling strength for each dimension; if numel=1 and dims > 1 then
    % we use that coupling strength for all dimensions
    coupling = [];
    
    % scale parameter; applied when prior is unspecified and in
    % initialization
    % maximimization of the model evidence using the training data
    % is performed to select a scale
    % parameter for the scales in maxevid under the assumption of a
    % 1./scale prior (uniform in the log domain)
    scale  = 1;
    
    % mask that can be used to access only a subset of a volume of data
    % when specifying the coupling; specify mask as size of input
    % dimensions
    mask = [];
    
    % some mutable options for laplacedegenerate_ep
    
    fraction    = 0.9;   % fraction or power for fractional/power EP
    niter       = 100;    % maximum number of iterations
    temperature = 1;      % forces MAP like behaviour for t->0
    tolerance   = 1e-5;   % convergence criterion
    degenerate  = [];     % whether or not to run in degenerate mode
    
    % learned parameters
    
    Gauss; % the EP estimate
    convergence; % whether or not EP converged
    logp; % approximate log model evidence
    
    % coupling strength for individual tasks; each feature is coupled to
    % its corresponding feature in the other tasks
    multitask = 0;
    ntasks = 1; % number of multitask components
    taskcoupling=100; % strong task coupling by default
    
    % mixed effects model (mixed=1 only couples fixed effects with
    % coupling; mixed=2 couples fixed and random effects with coupling)
    mixed = 0;
    nmixed = 0; % number of mixed effects components
    
  end
  
  methods
    
    function obj = blogreg(varargin)
      
      obj = obj@dml.method(varargin{:});
      
      if isempty(obj.precbias)
        % NOTE: EP doesn't seem to converge well for too large scales
        obj.precbias = max([1e2 obj.scale]);
      end
      
      if obj.taskcoupling > 0
        obj.taskcoupling = -obj.taskcoupling;
      end
      
      if isempty(obj.indims) && ~isempty(obj.mask)
        % infer input dimensions from optional mask
        obj.indims = size(obj.mask);
      end
      
    end
    
    function obj = train(obj,X,Y)
      
      % multiple outputs
      if ~iscell(Y) && size(Y,2) > 1
        obj = ft_mv_noutput('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % multitask learning
      if iscell(X) && (~obj.mixed || iscell(X{1}))
        obj.multitask = 1;
      end
      if obj.multitask, obj.ntasks = length(X); end
      
      % mixed effects
      if obj.mixed
        if obj.multitask
          obj.nmixed = length(X{1});
        else
          obj.nmixed = length(X);
        end
      end
      
      if ~obj.multitask && ~obj.mixed
        sz  = size(X);
        obj.nfeatures = prod(sz(2:end));
        nsamples = size(X,1);
      elseif (obj.multitask && ~obj.mixed) || (~obj.multitask && obj.mixed)
        sz  = size(X{1});
        obj.nfeatures = prod(sz(2:end));
        nsamples = mean(cellfun(@(x)(size(x,1)),X));
      elseif obj.multitask && obj.mixed
        sz  = size(X{1}{1});
        obj.nfeatures = prod(sz(2:end));
        nsamples = mean(cellfun(@(x)(size(x,1)),X{1})); % just an indication
      end
      
      if isempty(obj.degenerate)
        obj.degenerate = nsamples < obj.nfeatures;
      end
      if obj.verbose
        if obj.degenerate
          fprintf('running in degenerate mode\n');
        else
          fprintf('running in nondegenerate mode\n');
        end
      end
      
      % create prior; bias terms are added
      if isempty(obj.prior)
        obj.prior = obj.create_prior();
      elseif obj.verbose
        fprintf('using prespecified prior; bias term should be included in prior\n');
      end
      
      if obj.multitask || obj.mixed
        obj.prior = obj.transform_prior();
      end
      
      % change representation of the data
      [X,Y] = obj.transform_data(X,Y);
      
      % run the algorithm
      if isscalar(obj.scale)
        % learn model for fixed scale
        
        if obj.verbose
          fprintf('scaling prior with scale %g and bias scale %g\n',obj.scale,obj.precbias);
        end
        
        % scale prior local + bias term scaling
        obj.prior = obj.scale_prior(obj.scale);
        
        if obj.niter
          [obj.Gauss,terms,obj.logp,obj.convergence] = laplacedegenerate_ep(Y,X,obj.prior, ...
            'fraction',obj.fraction,'niter',obj.niter,'temperature',obj.temperature,'lambda',obj.scale,...
            'tol',obj.tolerance,'degenerate',obj.degenerate,'verbose',obj.verbose);
        end
        
      else
        
        lgp = -inf * ones(1,length(obj.scale));
        for j=1:length(obj.scale)
          
          if obj.verbose
            fprintf('scaling prior with scale %g and bias scale %g\n',obj.scale(j),obj.precbias);
          end
          tprior = obj.scale_prior(obj.scale(j));
          
          try
            
            if obj.niter
              [tgauss,terms,tlogp,tconvergence] = laplacedegenerate_ep(Y,X,tprior, ...
                'fraction',obj.fraction,'niter',obj.niter,'temperature',obj.temperature,'lambda',obj.scale(j),...
                'tol',obj.tolerance,'degenerate',obj.degenerate,'verbose',obj.verbose);
            end
            
            lgp(j) = tlogp - log(obj.scale(j)); % uniform prior on log scale
            if lgp(j) == max(lgp)
              
              gauss = tgauss;
              conv =  tconvergence;
              mprior = tprior;
              sc = obj.scale(j);
            end
            
          catch
            lgp(j) = -inf; % some error occurred
            disp(lasterr);
          end
          
        end
        
        obj.Gauss = gauss;
        obj.convergence = conv;
        obj.prior = mprior;
        obj.logp = lgp;
        obj.scale = sc;
        
        if obj.verbose
          fprintf('selected scale %f\n',sc);
        end
        
      end
      
      if obj.verbose
        if ~obj.convergence
          fprintf('EP did not converge\n');
        else
          fprintf('EP converged\n');
        end
      end
      
    end
    
    function Y = test(obj,X)
      
      G = obj.Gauss;
      
      D = obj.transform_data(X);
      
      % compute univariate means
      M = D * G.m;
      
      % compute univariate variances on the fly
      % also add the covariance for the betas to the output.
      
      nsamples = size(G.A,1);
      totsamples = size(D,1);
      
      scaledA = G.A .* (repmat(1./G.diagK',nsamples,1));
      W1 = G.A*scaledA';
      % now W1 = A * diag(1./obj.Gauss.diagK) * A'
      
      % add Delta (aka hatK)
      W1(1:(nsamples+1):numel(W1)) = W1(1:(nsamples+1):numel(W1)) + (1./G.hatK)';
      
      W2 = D*scaledA';
      % now W2 = D * diag(1./obj.Gauss.diagK) * A'
      
      scaledX = D .* (repmat(1./G.diagK',totsamples,1));
      W3 = D*scaledX';
      % now W3 = D * diag(1./obj.Gauss.diagK) * X'
      
      C = diag(W3 - (W2 / W1) * W2');
      % changed W2 * inv(W1) * W2' to (W2 / W1) * W2'
      
      if strcmp(obj.approximation,'probit') % probit approximation
        
        z = M .* (1 + pi .* C./8).^(-0.5);
        y = 1 ./ (1 + exp(-z));
        
      else % Gaussian quadrature to compute sigma(z) * N(z,m,c)
        
        nweights = obj.nweights;
        whermite = nan;
        
        % repeat until we have valid points
        while any(isnan(whermite(:)))
          
          [xhermite,whermite] = gausshermite(nweights);
          xhermite = xhermite(:);
          whermite = whermite(:);
          nhermite = length(whermite);
          
          nweights = nweights/2;
        end
        
        x = repmat(M,1,nhermite) + sqrt(C)*xhermite';
        g = logist(x);   % returns - log (1 + exp(-x)) with special attention for very small and very large x
        
        h = g + log(repmat(whermite',totsamples,1));
        maxh = max(h,[],2);
        
        y = exp(maxh) .* sum(exp(h - repmat(maxh,[1 size(h,2)])),2);
        
      end
      
      post = [y 1 - y];
      
      if ~obj.multitask && ~obj.mixed
        
        Y = post;
        
      elseif obj.multitask && ~obj.mixed
        
        Y = cell(1,obj.ntasks);
        nidx = 0;
        for j=1:obj.ntasks
          
          Y{j} = post((nidx+1):(nidx+size(X{j},1)),:);
          nidx = nidx + size(X{j},1);
          
        end
        
      elseif ~obj.multitask && obj.mixed
        
        Y = cell(1,obj.nmixed);
        nidx = 0;
        for j=1:obj.nmixed
          
          Y{j} = post((nidx+1):(nidx+size(X{j},1)),:);
          nidx = nidx + size(X{j},1);
          
        end
        
      elseif obj.multitask && obj.mixed
        
        Y = cell(1,obj.ntasks);
        nidx = 0;
        for j=1:obj.ntasks
          
          Y{j} = cell(1,obj.nmixed);
          for m=1:obj.nmixed
            
            Y{j}{m} = post((nidx+1):(nidx+size(X{j}{m},1)),:);
            nidx = nidx + size(X{j}{m},1);
            
          end
          
        end
        
      end
      
    end
    
    function m = model(obj)
      % MODEL returns the following parameters:
      %
      % m.importance importance values (posterior  - prior variance of auxiliary variables)
      % m.posterior posterior variance of auxiliary variables 
      % m.prior prior variance of auxiliary variables
      % m.bmean means of the regression coefficients; positive values indicate condition one
      % m.bvar variance of the betas
      %
      % In case of multiple tasks, each parameter is repeated for each task
      
      G = obj.Gauss;
      
      if ~(obj.multitask && obj.mixed)
        
        if ~obj.multitask && ~obj.mixed
          
          ntasks = obj.ntasks;
          N = (numel(G.auxC) - ntasks) / ntasks; % number of features
          offset = 0:(N+1):((ntasks-1)*(N+1));
          
        elseif obj.multitask && ~obj.mixed
          
          ntasks = obj.ntasks;
          N = (numel(G.auxC) - ntasks) / ntasks; % number of features
          offset = 0:(N+1):((ntasks-1)*(N+1));
          
        elseif ~obj.multitask && obj.mixed
          
          ntasks = obj.nmixed + 1;
          N = (numel(G.auxC) - ntasks + 1) / ntasks; % number of features
          offset = [0 N N + ((N+1):(N+1):((ntasks-2)*(N+1)))];
          
        end
        
        m = [];
        m.importance = cell(1,1+obj.nmixed);
        m.posterior = cell(1,1+obj.nmixed);
        m.prior = cell(1,1+obj.nmixed);
        m.bmean = cell(1,1+obj.nmixed);
        m.bvar = cell(1,1+obj.nmixed);
        
        for j=1:ntasks
          
          % Note: bias term must be included in prior
          varaux = G.auxC(offset(j) + (1:N)); % ignore bias term
          
          % compute variances of the auxiliary variables under the prior
          [L,dummy,S] = chol(sparse(obj.prior),'lower');
          invA = fastinvc(L);
          varprior = full(diag(S*invA*S'));
          varprior = varprior(offset(j) + (1:N));
          
          % mean and variance of the regression coefficients
          meanbeta = G.m(offset(j) + (1:N));
          varbeta = G.diagC(offset(j) + (1:N));
          
          % the mean and variance can be used to create eg 100(1-alpha)% credible
          % interval: ci(k,:) = [ ...
          % mu(k) - norminv(1-alpha/2)*sigma(k) ...
          % mu(k) + norminv(1-alpha/2)*sigma(k)
          % ]
          % with mu(k) = meanbeta(k) and sigma(k) = sqrt(varbeta(k)) such that
          % the alpha-importance map is defined as M_alpha = ci(:,1) < 0 & ci(:,2) > 0
          
          % model is posterior variance divided by prior variance of the
          % auxiliary variables; chose minus because of interpretation
          % problems...
          
          m.importance{j} = (varaux - varprior);
          m.posterior{j} = varaux;
          m.prior{j} = varprior;
          m.bmean{j} = meanbeta;
          m.bvar{j} = varbeta;
          
          if ~isempty(obj.mask)
            
            msk = obj.mask(:) ~= 0;
            mm = nan(obj.indims);
            
            mm(msk) = m.importance{j}; m.importance{j} = mm;
            mm(msk) = m.posterior{j}; m.posterior{j} = mm;
            mm(msk) = m.prior{j}; m.prior{j} = mm;
            mm(msk) = m.bmean{j}; m.bmean{j} = mm;
            mm(msk) = m.bvar{j}; m.bvar{j} = mm;
            
          end
          
        end
        
        % for mixed effects models, the first model is the fixed effects
        % model!
      
      else % obj.multitask && obj.mixed
        
        N1 = numel(G.auxC)./obj.ntasks;
        N = (N1-obj.nmixed) ./ (obj.nmixed+1);
        
        m = [];
        m.importance = cell(obj.ntasks,1+obj.nmixed);
        m.posterior = cell(obj.ntasks,1+obj.nmixed);
        m.prior = cell(obj.ntasks,1+obj.nmixed);
        m.bmean = cell(obj.ntasks,1+obj.nmixed);
        m.bvar = cell(obj.ntasks,1+obj.nmixed);
        
        for j=1:obj.ntasks
                    
          for k=1:(1+obj.nmixed)
            
            offset = (j-1)*N1 + (k-1)*(N+1);
            if k>1, offset = offset-1; end
            
            varaux = G.auxC(offset + (1:N)); % ignore bias term
            
            [L,dummy,S] = chol(sparse(obj.prior),'lower');
            invA = fastinvc(L);
            varprior = full(diag(S*invA*S'));
            varprior = varprior(offset + (1:N));
            meanbeta = G.m(offset + (1:N));
            varbeta = G.diagC(offset + (1:N));
            
            m.importance{j,k} = (varaux - varprior);
            m.posterior{j,k} = varaux;
            m.prior{j,k} = varprior;
            m.bmean{j,k} = meanbeta;
            m.bvar{j,k} = varbeta;
            
            if ~isempty(obj.mask)
              
              msk = obj.mask(:) ~= 0;
              mm = nan(obj.indims);
              
              mm(msk) = m.importance{j,k}; m.importance{j,k} = mm;
              mm(msk) = m.posterior{j,k}; m.posterior{j,k} = mm;
              mm(msk) = m.prior{j,k}; m.prior{j,k} = mm;
              mm(msk) = m.bmean{j,k}; m.bmean{j,k} = mm;
              mm(msk) = m.bvar{j,k}; m.bvar{j,k} = mm;
                          
            end
            
          end
          
        end
        
      end
      
      if numel(m.importance)==1
        m.importance = m.importance{1};
        m.posterior = m.posterior{1};
        m.prior = m.prior{1};
        m.bmean = m.bmean{1};
        m.bvar = m.bvar{1};
      end
      
    end
    
    function [Y,X,B,P,u,v] = sample(obj,M)
      % This function samples from the prior and creates M vectors for
      % the regression coefficients beta; this function replicates
      % sample_betas and then computes a dataset from it.
      %
      % Y   : outputs
      % X   : covariates
      % B   : regression coefficients
      % P   : Bernoulli probabilities
      % u,v : auxiliary variables
      %
      % e.g.,
      %
      % b = blogreg('scale',0.1)
      % [Y,X,B,P,u,v] = sample(b,10000);
      % a1 = histc(B(1,:),linspace(-1,1,20));
      
      if obj.verbose
        fprintf('sampling betas from auxiliary variables using scaled prior\n');
      end
      
      if ~isfield(obj,'prior')
        % create univariate prior with a given scale
        obj.prior = obj.create_prior();
        pri = obj.scale_prior(obj.scale);
      else
        pri = obj.prior;
      end
      
      if nargin < 2, M = 1; end
      n = size(pri,1);
      
      % get samples for auxiliary variables
      u = sample_from_prior(zeros(n,1),pri,M);
      v = sample_from_prior(zeros(n,1),pri,M);
      
      % get samples for betas
      B = normrnd(zeros(n,M),sqrt(u.^2 + v.^2));
      
      % create random dataset
      X = randn(size(B));
      X(size(B,1),:) = 1; % bias
      
      P = 1 ./ (1 + exp(-sum(X .* B)));
      
      Y = (P < rand(size(P))) + 1;
      
      % data as nexamples X nfeatures
      X = X(1:(size(B,1)-1),:)'; % ignore bias term
      Y = Y';
      
    end
    
    function [tdata,tdesign] = transform_data(obj,X,Y)
      % recreate dataset and design matrix
      
      if ~obj.multitask && ~obj.mixed
        
        % add bias term
        tdata = [X(1:size(X,1),:) ones(size(X,1),1)];
        
        % transform design to +1/-1 representation
        if nargin==3, tdesign = 3-2*Y; end
        
      elseif obj.multitask && ~obj.mixed
        
        if nargin==3,
          tdesign = cell2mat(Y(:));
          tdesign = 3-2*tdesign(:,1);
        end
        
        % add bias term
        for j=1:obj.ntasks, X{j} = [X{j}(1:size(X{j},1),:) ones(size(X{j},1),1)]; end
        
        tdata = blkdiag(X{:});
        
      elseif ~obj.multitask && obj.mixed
        
        if nargin==3,
          tdesign = cell2mat(Y(:));
          tdesign = 3-2*tdesign(:,1);
        end
        
        % add bias term
        for j=1:obj.nmixed, X{j} = [X{j}(1:size(X{j},1),:) ones(size(X{j},1),1)]; end
        
        % fixed effects without bias
        fixed = cell2mat(X(:));
        fixed = fixed(:,1:(end-1));
        
        tdata = [fixed blkdiag(X{:})];
        
      else % obj.multitask && obj.mixed
        
        if nargin==3,
          tdesign = [];
          for k=1:length(Y)
            td = cell2mat(Y{k}(:));
            td = 3-2*td(:,1);
            tdesign = cat(1,tdesign,td);
          end
        end
        
        tdata = cell(1,length(X));
        for k=1:length(X)
          
          % add bias term
          XX = cell(size(X));
          for j=1:obj.nmixed, XX{j} = [X{k}{j}(1:size(X{k}{j},1),:) ones(size(X{k}{j},1),1)]; end
          
          % fixed effects without bias
          fixed = cell2mat(XX(:));
          fixed = fixed(:,1:(end-1));
          
          tdata{k} = [fixed blkdiag(XX{:})];
          
        end
        
        tdata = blkdiag(tdata{:});
        
      end
      
      tdata = double(tdata);
      if nargin==3, tdesign = double(tdesign); end
      
    end
    
  end
  
  methods(Access=protected)
    
    
    function prior = create_prior(obj)
      
      nf = obj.nfeatures;
      
      if isempty(obj.coupling) || ~any(obj.coupling)
        
        if obj.verbose
          fprintf('using decoupled prior\n');
        end
        
        prior = spalloc(nf,nf,nf+1);
        prior(1:(nf+1):numel(prior)) = 1;
        
        if ~isempty(obj.mask)
          prior = prior(find(obj.mask(:)),find(obj.mask(:)));
        end
        
      else
        
        assert(~isempty(obj.indims));
        
        if numel(obj.coupling)==1 && numel(obj.indims) > 1
          obj.coupling = obj.coupling*ones(1,numel(obj.indims));
        end
        
        obj.coupling(obj.coupling > 0) = -obj.coupling(obj.coupling > 0);
        
        if obj.verbose
          fprintf('using prior with coupling [');
          for j=1:length(obj.coupling)
            fprintf(' %g',obj.coupling(j));
          end
          fprintf(' ]\n');
        end
        
        prior = construct_prior(obj.indims,obj.coupling,'mask',find(obj.mask(:)),'circulant',[0 0 0 0]);
        
      end
      
      prior(nf+1,nf+1) = 1;
      
    end
    
    function prior = transform_prior(obj)
      % redefine prior for multitask and/or mixed effects models
      
      nf = obj.nfeatures;
      
      % recreate prior
      prior = obj.prior;
      
      if obj.multitask && ~obj.mixed
        
        pp = repmat({prior},[1 obj.ntasks]);
        prior = blkdiag(pp{:});
        
        % add task coupling for multitask learning
        if obj.multitask
          
          blk = obj.ntasks*(nf+1)^2;
          didx = (1:(size(prior,1)+1):(size(prior,1)*(nf)));
          cc = repmat(obj.taskcoupling,[1 length(didx)]);
          for j=1:obj.ntasks
            for k=(j+1):obj.ntasks
              prior(((j-1)*blk+(k-1)*(nf+1))+didx) = cc;
              prior(((k-1)*blk+(j-1)*(nf+1))+didx) = cc;
            end
          end
        end
        
      elseif ~obj.multitask && obj.mixed
        
        if obj.mixed==2 % couple fixed and random effects
          
          P = prior;
          pp = cat(2,{prior(1:nf,1:nf)},repmat({P},[1 obj.nmixed]));
          prior = blkdiag(pp{:});
          
        else % couple fixed effects only
          
          P = zeros(nf+1,nf+1);
          P(1:(nf+2):end) = prior(1:(nf+2):end);
          
          pp = cat(2,{prior(1:nf,1:nf)},repmat({P},[1 obj.nmixed]));
          prior = blkdiag(pp{:});
          
        end
        
      else % obj.multitask && obj.mixed
        
        if obj.mixed==2 % couple fixed and random effects
          
          P = prior;
          pp = cat(2,{prior(1:nf,1:nf)},repmat({P},[1 obj.nmixed]));
          prior = blkdiag(pp{:});
          
        else % couple fixed effects only
          
          P = zeros(nf+1,nf+1);
          P(1:(nf+2):end) = prior(1:(nf+2):end);
          
          pp = cat(2,{prior(1:nf,1:nf)},repmat({P},[1 obj.nmixed]));
          prior = blkdiag(pp{:});
          
        end
        
        pp = repmat({prior},[1 obj.ntasks]);
        prior = blkdiag(pp{:});
        
        % add task coupling for multitask learning
        if obj.multitask
          
          nf2 = (nf+1)*(obj.nmixed+1)-1;
          
          blk = obj.ntasks*nf2^2;
          if obj.mixed==1 % couple fixed effects only
            didx = (1:(size(prior,1)+1):(size(prior,1)*nf));
          else % couple fixed and random effects
            didx = (1:(size(prior,1)+1):(size(prior,1)*nf2));
          end
          cc = repmat(obj.taskcoupling,[1 length(didx)]);
          for j=1:obj.ntasks
            for k=(j+1):obj.ntasks
              prior(((j-1)*blk+(k-1)*nf2)+didx) = cc;
              prior(((k-1)*blk+(j-1)*nf2)+didx) = cc;
            end
          end
        end
        
        % remove coupling of bias term
        for j=1:obj.ntasks
          for k=(j+1):obj.ntasks
            for m=1:obj.nmixed
              prior((j-1)*nf2 + nf + m*(nf+1),(k-1)*nf2 + nf + m*(nf+1)) = 0;
              prior((k-1)*nf2 + nf + m*(nf+1),(j-1)*nf2 + nf + m*(nf+1)) = 0;
            end
          end
        end
        
      end
      
    end
    
    function KK = scale_prior(obj,scale)
      % SCALE_PRIOR scales a prior precision matrix
      
      nfeatures = size(obj.prior,1);
      
      obj.prior(1:(nfeatures+1):numel(obj.prior)) = 0;
      obj.prior(1:(nfeatures+1):numel(obj.prior)) = 1 - sum(obj.prior,2);
      
      [L,d1,d2] = chol(obj.prior,'lower');
      
      if d1
        error('matrix not p.d.!');
      else
        C = speye(nfeatures) .* fastinvc(L); % only get diagonal elements
        C     = d2*C*d2'; % reorder
      end
      
      Csqrt = sqrt(C);
      
      KK = (Csqrt*obj.prior*Csqrt)./sqrt(scale(:)*scale(:)'); % suitable for mixed effects
      
      % add bias term scaling
      if ~obj.multitask && ~obj.mixed
        
        KK(end,end) = 1/obj.precbias;
        
      elseif obj.multitask && ~obj.mixed
        
        for j=1:obj.ntasks
          KK(j*(obj.nfeatures+1),j*(obj.nfeatures+1)) = 1/obj.precbias;
        end
        
      elseif ~obj.multitask && obj.mixed
        
        % no bias for fixed effects
        for j=1:obj.ntasks
          KK(obj.nfeatures + j*(obj.nfeatures+1),obj.nfeatures + j*(obj.nfeatures+1)) = 1/obj.precbias;
        end
        
      elseif obj.multitask && obj.mixed
        
        nf = obj.nfeatures;
        nf2 = (nf+1)*(obj.nmixed+1)-1;
        
        for j=1:obj.ntasks
          for m=1:obj.nmixed
            KK((j-1)*nf2 + nf + m*(nf+1),(j-1)*nf2 + nf + m*(nf+1)) = 1/obj.precbias;
          end
        end
        
      end
      
    end
    
  end
  
end
