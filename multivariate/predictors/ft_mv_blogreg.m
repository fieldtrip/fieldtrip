classdef ft_mv_blogreg < ft_mv_predictor
%FT_MV_BLOGREG Bayesian logistic regression with spatiotemporal interactions and
%possibility for transfer learning
%
% The scale property specifies the regularization. The bigger the scale,
% the less the regression coefficients will be regularized towards zero.
%
% ft_mv_blogreg allows feature to become coupled. This is done through the
% coupling property, which specifies for each input dimension the coupling
% strength in that dimension. Also, the input dimensions of the original
% data must be specified through the indims property.
% In case the input data X is just a region of interest then we can use a
% mask to indicate which features from the original volume of data are
% represented by X. This still allows us to use the above approach to
% specify the coupling.
%
% Transfer learning is implemented by augmenting the data matrix as
% [ subject1data      0       ]
% [      0       subject2data ]
% etc.
%
% NOTE: 
%   a bias term is added to the model and should not be included explicitly
%
% Refs:
% van Gerven et al. Efficient Bayesian multivariate fMRI analysis using a sparsifying
% spatio-temporal prior. Neuroimage, 2010
% van Gerven & Simanova. Concept classification with Bayesian multitask
% learning, NAACL, 2010
%
% Copyright (c) 2010, Marcel van Gerven
%

    properties

      % number of features
      nfeatures
      
      % coupling strength for individual tasks; each feature is coupled to
      % its corresponding feature in the other tasks
      taskcoupling=100; % strong coupling by default
      ntasks = 1; % no transfer learning by default

      % input dimensions of the data; is used to generate prior
      indims = [];
     
       % prior precision matrix of the auxiliary variables
      prior
            
      % scale of the bias term (bias term will be appended to the model)
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
     
    end

    methods
       
       function obj = ft_mv_blogreg(varargin)
      
        obj = obj@ft_mv_predictor(varargin{:});
        
        if isempty(obj.precbias)
          % NOTE: EP doesn't seem to converge well for too large scales
          obj.precbias = max([1e2 obj.scale]);
        end
        
        if obj.taskcoupling > 0
          obj.taskcoupling = -obj.taskcoupling;
        end
        
       end
       
       function obj = train(obj,X,Y)
         
         % missing data
         if any(isnan(X(:))) || any(isnan(Y(:))), error('method does not handle missing data'); end
         
         % multiple outputs
         if size(Y,2) > 1
           obj = ft_mv_noutput('mvmethod',obj);
           obj = obj.train(X,Y);
           return;
         end
         
         % transfer learning
         transfer = iscell(X);
         
         if transfer
           obj.ntasks = length(X);
           obj.nfeatures = size(X{1},2);
           nsamples = size(X{1},1);
          
         else
           obj.nfeatures = size(X,2);
           nsamples = size(X,1);
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
         
         if isempty(obj.prior)
           obj.prior = obj.create_prior();
         elseif obj.verbose
           fprintf('using prespecified prior\n');
         end
         
         if transfer % transfer learning
           obj.prior = obj.transform_prior();
         end
         
         [X,Y] = obj.transform_data(X,Y);
           
         % run the algorithm
         if isscalar(obj.scale)
           % learn model for fixed scale
           
           if obj.verbose
             fprintf('scaling prior with scale %g and bias scale %g\n',obj.scale,obj.precbias);
           end
           
           obj.prior = scale_prior(obj.prior,'lambda',obj.scale);
           
           if transfer
             for j=1:obj.ntasks
               obj.prior(j*(obj.nfeatures+1),j*(obj.nfeatures+1)) = obj.precbias;
             end
           elseif size(obj.prior,1)~=size(X,2)
             % add bias term if not yet done
             obj.prior(size(X,2),size(X,2)) = obj.precbias;
           end
           
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
             tprior = scale_prior(obj.prior,'lambda',obj.scale(j));
             
             if transfer
               for j=1:obj.ntasks
                 tprior(j*(obj.nfeatures+1),j*(obj.nfeatures+1)) = obj.precbias;
               end
             elseif size(tprior,1)~=size(X,2)
               % add bias term if not yet done
               tprior(size(X,2),size(X,2)) = obj.precbias;
             end
 
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
         if obj.ntasks > 1
           totsamples = sum(cellfun(@(x)(size(x,1)),X));
         else
           totsamples = size(X,1);
         end
         
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
         
         if obj.ntasks == 1
           Y = post;
         else
           Y = cell(1,obj.ntasks);
           nidx = 0;
           for j=1:obj.ntasks
             
             Y{j} = post((nidx+1):(nidx+size(X{j},1)),:);
             nidx = nidx + size(X{j},1);
             
           end
         end
         
       end
       
       function [tdata,tdesign] = transform_data(obj,X,Y)
         % recreate dataset and design matrix
         
         if obj.ntasks == 1
           
           % add bias term
           tdata = [X ones(size(X,1),1)];
           
           % transform design to +1/-1 representation
           if nargin==3, tdesign = 3-2*Y; end
           
         else % transfer learning
           
           nf = obj.nfeatures;
           totsamples = sum(cellfun(@(x)(size(x,1)),X));
           
           % sparse might be faster for many tasks...
           tdata = zeros(totsamples,obj.ntasks*(nf+1));
           if nargin==3, tdesign = zeros(totsamples,1); end
           
           nidx = 0; fidx = 0;
           for j=1:obj.ntasks
             
             tdata((nidx+1):(nidx+size(X{j},1)),(fidx+1):(fidx+nf+1)) = [X{j} ones(size(X{j},1),1)];
             if nargin==3
               tdesign((nidx+1):(nidx+size(X{j},1))) = 3-2*Y{j}(:,1);
             end
             
             nidx = nidx + size(X{j},1);
             fidx = fidx + nf + 1;
             
           end
           
         end
         
       end
       
       function prior = transform_prior(obj)
         % recreate prior for transfer learning
         
         nf = obj.nfeatures;
         
         % recreate prior
         prior = obj.prior;
         
         if size(prior,1)==nf
           % add bias term if not yet done
           prior(nf+1,nf+1) = 1;
         end
         
         pp = repmat({prior},[1 obj.ntasks]);
         prior = blkdiag(pp{:});
         
         % add task coupling
         if obj.taskcoupling
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
         
       end
       
       function [m,desc] = model(obj)
         % return the variances of the auxiliary variables as the model; this
         % determines in turn the magnitude of the betas through: U = u^2 + v^2
         % we output variances relative to the prior variances
         %
         % other return values:
         % mean of the betas
         % variance of the betas
         %
         
         % Note: bias term must be included in prior
         varaux = obj.Gauss.auxC(1:(end-1)); % ignore bias term
         
         % compute variances of the auxiliary variables under the prior
         [L,dummy,S] = chol(sparse(obj.prior),'lower');
         invA = fastinvc(L);
         varprior = full(diag(S*invA*S'));
         varprior = varprior(1:(end-1));
         
         % mean and variance of the regression coefficients
         meanbeta = obj.Gauss.m(1:(end-1));
         varbeta = obj.Gauss.diagC(1:(end-1));
         
         % model is posterior variance divided by prior variance of the
         % auxiliary variables; chose minus because of interpretation
         % problems...
         
         m = cell(5,1);
         m{1} = (varaux - varprior);
         m{2} = varaux;
         m{3} = varprior;
         m{4} = meanbeta;
         m{5} = varbeta;
         
         if ~isempty(obj.mask)
           
           msk = obj.mask(:) ~= 0;
           mm = nan(obj.indims);
           for j=1:5
             mm(msk) = m{j};
             m{j} = mm;
           end
           
         end
         
         desc = { ...
           'importance values (posterior  - prior variance of auxiliary variables)' ...
           'posterior variance of auxiliary variables' ...
           'prior variance of auxiliary variables' ...
           'means of the regression coefficients; positive values indicate condition one' ...
           'variance of the regression coefficients' }';
         
         % split for each task
         if obj.ntasks > 1
           mm = cell(5,obj.ntasks);
           for c=1:5
             
             % re-add bias term removed by blogreg
             m{c} = [m{c}; 1];
             
             tmp = reshape(m{c},[numel(m{c})./obj.ntasks obj.ntasks]);
             
             % remove bias
             tmp = tmp(1:(end-1),:);
             
             for j=1:obj.ntasks
               mm{c,j} = tmp(:,j);
             end
           end
           m = mm;
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
         
         if ~isfield(obj.params,'prior')
           % create univariate prior with a given scale
           pri = obj.create_prior(1);
           pri = scale_prior(pri,'lambda',obj.scale);
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
       
       
    end
    
    methods(Access=protected)
      
      function prior = create_prior(obj)
        
        nf = obj.nfeatures;
        
        if isempty(obj.coupling) || isempty(obj.indims) || ~any(obj.coupling)
          
          if obj.verbose
            fprintf('using decoupled prior\n');
          end
          
          prior = spalloc(nf,nf,nf+1);
          prior(1:(nf+1):numel(prior)) = 1;
          
          if ~isempty(obj.mask)
            prior = prior(find(obj.mask(:)),find(obj.mask(:)));
          end
          
        else
          
          if numel(obj.coupling)==1 && numel(obj.indims) > 1
            obj.coupling = obj.coupling*ones(1,numel(obj.indims));
          end
          
          obj.coupling(obj.coupling > 0) = -obj.coupling(obj.coupling > 0);
          
          if obj.verbose
            fprintf('using prior with coupling [ ');
            fprintf('%g]\n',obj.coupling);
          end
          
          prior = construct_prior(obj.indims,obj.coupling,'mask',find(obj.mask(:)),'circulant',[0 0 0 0]);
          
        end
        
      end
      
    end
    
end
