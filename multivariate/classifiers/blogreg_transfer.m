classdef blogreg_transfer < blogreg & transfer_learner
%Bayesian logistic regression with spatiotemporal interactions;
% transfer learning formulation. 
%
% Transfer learning can be implemented for any classifier that couples
% variables by augmenting the data matrix as
%
% [ subject1data      0       ]
% [      0       subject2data ]
%
% etc.
%
% assuming nfeatures is the same for all data; we might change
% this in the future.
%         
% PARAMETERS:
% p.prior; % the used prior
% p.Gauss; % the EP estimate
% p.convergence; % whether or not EP converged
% p.logp; % approximate log model evidence
% p.scale; % selected scale (in case of multiple scales); large scale is strong regularization!
% p.ntasks; % number of used tasks
%
% EXAMPLE:
% [a,b,c] = test_procedure({standardizer blogreg_transfer('taskcoupling',100)},0.8)
% bar([c.model{1,1} c.model{1,2}])
%
% Copyright (c) 2010, Marcel van Gerven
%

    properties

      % coupling strength for individual tasks; each feature is coupled to
      % its corresponding feature in the other tasks
      taskcoupling=100; % strong coupling by default
      
    end

    methods
       
      function obj = blogreg_transfer(varargin)
      
           obj = obj@blogreg(varargin{:});
                 
           if obj.taskcoupling > 0
             obj.taskcoupling = -obj.taskcoupling;
           end
           
       end
       
       function p = estimate(obj,X,Y)
                  
         % take out missing data
         for c=1:length(X)
           lab = mvmethod.labeled(Y{c});
           X{c} = X{c}(lab,:);
           Y{c} = Y{c}(lab,:);
         end
         
         if iscell(obj.indims)
           obj.dims = obj.indims{1};
         else
           obj.dims = obj.indims;
         end

         nfeatures = size(X{1},2);
         
         p.ntasks = length(X);
         
         if isempty(obj.prior)
           p.prior = obj.create_prior(nfeatures);
         else           
           if obj.verbose
             fprintf('using prespecified prior\n');
           end
           p.prior = obj.prior;
         end
           
         if isempty(obj.degenerate)
           obj.degenerate = size(X{1},1) < nfeatures;
         end
         
         if obj.verbose
           if obj.degenerate
             fprintf('running in degenerate mode\n');
           else
             fprintf('running in nondegenerate mode\n');
           end
         end
                 
         % recreate dataset and design matrix
         
         totsamples = sum(cellfun(@(x)(size(x,1)),X));
         
         % sparse might be faster for many tasks...
         tdata = zeros(totsamples,p.ntasks*(nfeatures+1));
         tdesign = zeros(totsamples,1);
         
         nidx = 0; fidx = 0;
         for j=1:p.ntasks
           
           tdata((nidx+1):(nidx+size(X{j},1)),(fidx+1):(fidx+nfeatures+1)) = [X{j} ones(size(X{j},1),1)];
           
           tdesign((nidx+1):(nidx+size(X{j},1))) = 3-2*Y{j}(:,1);
           
           nidx = nidx + size(X{j},1);
           fidx = fidx + nfeatures + 1;
           
         end
         
         clear X;
         clear Y;
         
         % recreate prior
         
         if size(p.prior,1)==nfeatures
           % add bias term if not yet done
           p.prior(nfeatures+1,nfeatures+1) = 1;
         end
         
         pp = repmat({p.prior},[1 p.ntasks]);
         p.prior = blkdiag(pp{:});
         
         % add task coupling
         if obj.taskcoupling
           blk = p.ntasks*(nfeatures+1)^2;
           didx = (1:(size(p.prior,1)+1):(size(p.prior,1)*(nfeatures)));
           cc = repmat(obj.taskcoupling,[1 length(didx)]);
           for j=1:p.ntasks
             for k=(j+1):p.ntasks
               p.prior(((j-1)*blk+(k-1)*(nfeatures+1))+didx) = cc;
               p.prior(((k-1)*blk+(j-1)*(nfeatures+1))+didx) = cc;
             end
           end
         end
         
         % run the algorithm
          
         if isscalar(obj.scale)
           % learn model for fixed scale
           
           p.prior = scale_prior(p.prior,'lambda',obj.scale);           
           
           % override for the bias terms
           for j=1:p.ntasks
             p.prior(j*(nfeatures+1),j*(nfeatures+1)) = obj.precbias;
           end
           
           if obj.niter
             [p.Gauss,terms,p.logp,p.convergence] = laplacedegenerate_ep(tdesign,tdata,p.prior, ...
               'fraction',obj.fraction,'niter',obj.niter,'temperature',obj.temperature,'lambda',obj.scale,...
               'tol',obj.tolerance,'degenerate',obj.degenerate,'verbose',obj.verbose);
           end
           
         else
           
           lgp = -inf * ones(1,length(obj.scale));
           for j=1:length(obj.scale)
             
             tprior = scale_prior(p.prior,'lambda',obj.scale(j));
             
             % override for the bias terms
             for jj=1:p.ntasks
               tprior(jj*(nfeatures+1),jj*nfeatures+1) = obj.precbias;
             end
             
             try
               
               if obj.niter
                 [p.Gauss,terms,p.logp,p.convergence] = laplacedegenerate_ep(tdesign,tdata,tprior, ...
                   'fraction',obj.fraction,'niter',obj.niter,'temperature',obj.temperature,'lambda',obj.scale(j),...
                   'tol',obj.tolerance,'degenerate',obj.degenerate,'verbose',obj.verbose);
               end
               
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
       
       function Y = map(obj,X)       
         
         G = obj.params.Gauss;
         
         % recreate dataset and design matrix
         
         nfeatures = size(X{1},2);
         
         totsamples = sum(cellfun(@(x)(size(x,1)),X));
         
         D = zeros(totsamples,obj.params.ntasks*(nfeatures+1));
         
         nidx = 0; fidx = 0;
         for j=1:obj.params.ntasks
           
           D((nidx+1):(nidx+size(X{j},1)),(fidx+1):(fidx+nfeatures+1)) = [X{j} ones(size(X{j},1),1)];
           
           nidx = nidx + size(X{j},1);
           fidx = fidx + nfeatures + 1;
           
         end
         
         % compute univariate means
         M = D * G.m;
         
         % compute univariate variances on the fly
         % also add the covariance for the betas to the output.
         
         nsamples = size(G.A,1);
         
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
         
         Y = cell(1,obj.params.ntasks);
         nidx = 0;
         for j=1:obj.params.ntasks
           
           Y{j} = post((nidx+1):(nidx+size(X{j},1)),:);           
           nidx = nidx + size(X{j},1);
           
         end
         
         
       end
       
       function [m,desc] = getmodel(obj)
         
         [mm,desc] = obj.getmodel@blogreg();

         % split for each task
         m = cell(5,obj.params.ntasks);
         for c=1:5

           % re-add bias term removed by blogreg
           mm{c} = [mm{c}; 1];
           
           tmp = reshape(mm{c},[numel(mm{c})./obj.params.ntasks obj.params.ntasks]);
           
           % remove bias 
           tmp = tmp(1:(end-1),:);
           
           for j=1:obj.params.ntasks
             m{c,j} = tmp(:,j);
           end
         end
         
       end
       
       
    end
    
end
