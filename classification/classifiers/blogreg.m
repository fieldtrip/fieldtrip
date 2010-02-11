classdef blogreg < classifier
%Bayesian logistic regression with spatiotemporal interactions
%
%
% PARAMETERS:
% p.prior; % the used prior
% p.Gauss; % the EP estimate
% p.convergence; % whether or not EP converged
% p.logp; % approximate log model evidence
% p.scale; % selected scale (in case of multiple scales)
%
% NOTE: 
%   a bias term is added to the model
%
% Copyright (c) 2009, Marcel van Gerven


    properties

      % prior precision matrix of the auxiliary variables
      prior
            
      % precision of the bias term (will be added to the model)
      precbias = []; % small precision translates into high variance
      
      % (either 'probit' or gaussian 'quadrature' to approximate posterior) 
      approximation = 'quadrature'      
      
      % number of weights for gaussian quadrature
      nweights = 100;
      
      % dimensions of the data; can be used to generate prior
      dims = [];
      
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
      % when specifying the coupling
      mask = [];
            
      % some mutable options for laplacedegenerate_ep
      
      fraction    = 0.9;   % fraction or power for fractional/power EP
      niter       = 100;    % maximum number of iterations
      temperature = 1;      % forces MAP like behaviour for t->0
      tolerance   = 1e-5;   % convergence criterion
      degenerate  = [];     % whether or not to run in degenerate mode
      
    end

    methods
       
      function obj = blogreg(varargin)
      
           % check availability
           if ~exist('laplacedegenerate_ep','file')
               error('not yet available...');
           end

           obj = obj@classifier(varargin{:});
      
           if isempty(obj.precbias)
                              
               % NOTE: EP doesnt seem to converge for low precisions 
               % (i.e., large scales)
               
               % we choose the precision of the bias term to  equal the scale
               obj.precbias = min(1e-2,1./obj.scale(1));
           end                                                                
           
       end
       
       function p = estimate(obj,X,Y)
                  
         p = obj.params;
         
         obj.dims = obj.indims(2:end);
                        
         if isempty(obj.prior)
           p.prior = obj.create_prior(size(X,2));
         else           
           if obj.verbose
             fprintf('using prespecified prior\n');
           end      
           p.prior = obj.prior;
         end
           
         if isempty(obj.degenerate)
           obj.degenerate = size(X,1) < size(X,2);
         end         
         if obj.verbose
           if obj.degenerate
             fprintf('running in degenerate mode\n');
           else
             fprintf('running in nondegenerate mode\n');
           end
         end
                  
         % add bias term
         X = [X ones(size(X,1),1)];
         
         % run the algorithm
         
         % transform design to +1/-1 representation
         design = 3-2*Y(:,1);
         
         if isscalar(obj.scale)
           % learn model for fixed scale
           
           p.prior = scale_prior(p.prior,'lambda',obj.scale);
           
           if size(p.prior,1)~=size(X,2)
             % add bias term if not yet done
             p.prior(size(X,2),size(X,2)) = obj.precbias;
           end
           
           [p.Gauss,terms,p.logp,p.convergence] = laplacedegenerate_ep(design,X,p.prior, ...
             'fraction',obj.fraction,'niter',obj.niter,'temperature',obj.temperature,'lambda',obj.scale,...
             'tol',obj.tolerance,'degenerate',obj.degenerate,'verbose',false);
           
          else
           
           lgp = -inf * ones(1,length(obj.scale));
           for j=1:length(obj.scale)
             
             tprior = scale_prior(p.prior,'lambda',obj.scale(j));
             
             if size(tprior,1)~=size(X,2)
               % add bias term if not yet done
               tprior(size(X,2),size(X,2)) = obj.precbias;
             end
             
             try 
               [p.Gauss,terms,p.logp,p.convergence] = laplacedegenerate_ep(design,X,tprior, ...
                 'fraction',obj.fraction,'niter',obj.niter,'temperature',obj.temperature,'lambda',obj.scale(j),...
                 'tol',obj.tolerance,'degenerate',obj.degenerate,'verbose',false);
               
               lgp(j) = p.logp - log(obj.scale(j)); % uniform prior on log scale
               
               if lgp(j) == max(lgp)
                 
                 gauss = p.Gauss;
                 conv =  p.convergence;
                 mprior = tprior;
                 sc = obj.scale(j);
                 
               end
               
             catch
              lgp(j) = -inf; % some error occurred
             end
           end
           
           p.Gauss = gauss;
           p.convergence = conv;
           p.prior = mprior;
           p.logp = lgp;
           p.scale = sc;
           
           if obj.verbose
             fprintf('selected scale %f\n',sc);
           end
           
         end
                           
         if obj.verbose
           if ~p.convergence
             fprintf('EP did not converge\n');
           else
             fprintf('EP converged\n');
           end
         end
          
       end
       
       function post = map(obj,X)       
         
         X = [X ones(size(X,1),1)];
         
         G = obj.params.Gauss;
         
         % compute univariate means
         M = X * G.m;
         
         % compute univariate variances on the fly
         % also add the covariance for the betas to the output.
         
         nsamples = size(G.A,1);
         
         scaledA = G.A .* (repmat(1./G.diagK',nsamples,1));
         W1 = G.A*scaledA';
         % now W1 = A * diag(1./obj.Gauss.diagK) * A'
         
         % add Delta (aka hatK)
         W1(1:(nsamples+1):numel(W1)) = W1(1:(nsamples+1):numel(W1)) + (1./G.hatK)';
         
         W2 = X*scaledA';
         % now W2 = X * diag(1./obj.Gauss.diagK) * A'
         
         scaledX = X .* (repmat(1./G.diagK',size(X,1),1));
         W3 = X*scaledX';
         % now W3 = X * diag(1./obj.Gauss.diagK) * X'
         
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
           
           h = g + log(repmat(whermite',size(X,1),1));
           maxh = max(h,[],2);
           
           y = exp(maxh) .* sum(exp(h - repmat(maxh,[1 size(h,2)])),2);
           
         end
         
         post = [y 1 - y];
         
       end
       
       function [m,desc] = getmodel(obj)
         % return the variances of the auxiliary variables as the model; this
         % determines in turn the magnitude of the betas through: U = u^2 + v^2
         % we output variances relative to the prior variances
         %
         % other return values:
         % mean of the betas
         % variance of the betas
         %
         
         % Note: bias term must be included in prior
         varaux = obj.params.Gauss.auxC(1:(end-1)); % ignore bias term
         
         % compute variances of the auxiliary variables under the prior
         [L,dummy,S] = chol(sparse(obj.params.prior),'lower');
         invA = fastinvc(L);
         varprior = full(diag(S*invA*S'));
         varprior = varprior(1:(end-1));
         
         % mean and variance of the regression coefficients
         meanbeta = obj.params.Gauss.m(1:(end-1));
         varbeta = obj.params.Gauss.diagC(1:(end-1));
         
         % model is posterior variance divided by prior variance of the
         % auxiliary variables; chose minus because of interpretation
         % problems...
                
         m = cell(5,1);
         m{1} = (varaux - varprior);
         m{2} = varaux;
         m{3} = varprior;
         m{4} = meanbeta;
         m{5} = varbeta; 
         
         desc = { ...
           'importance values (posterior  - prior variance of auxiliary variables)' ...
           'posterior variance of auxiliary variables' ...
           'prior variance of auxiliary variables' ...
           'means of the regression coefficients' ...
           'variance of the regression coefficients' }';
         
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
         
         if obj.verbose
           fprintf('sampling betas from auxiliary variables using scaled prior\n');
         end
         
         if nargin < 2, M = 1; end
         n = size(obj.params.prior,1);
         
         % get samples for auxiliary variables
         u = sample_from_prior(zeros(n,1),obj.params.prior,M);
         v = sample_from_prior(zeros(n,1),obj.params.prior,M);
         
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
      
      function prior = create_prior(obj,nfeatures)
        
        if isempty(obj.coupling) || isempty(obj.dims) || ~any(obj.coupling)
          
          if obj.verbose
            fprintf('using decoupled prior\n');
          end
          
          prior = spalloc(nfeatures,nfeatures,nfeatures+1);
          prior(1:(nfeatures+1):numel(prior)) = 1;
          
        else
          
          if numel(obj.coupling)==1 && numel(obj.dims) > 1
            obj.coupling = obj.coupling*ones(1,numel(obj.dims));
          end
          
          obj.coupling(obj.coupling > 0) = -obj.coupling(obj.coupling > 0);
          
          if obj.verbose
            fprintf('using prior with coupling [ ');
            fprintf('%g ',obj.coupling);
            fprintf('], scale %g, and bias precision %g\n',obj.scale,obj.precbias);
          end
          
          prior = construct_prior(obj.dims,obj.coupling,'mask',obj.mask,'circulant',[0 0 0 0]);
          
        end
        
      end
      
    end
    
end
