classdef gnb < classifier
%GNB generalized naive Bayes classifier
%
% If no distribution is specified then we assume a gaussian distribution
% and invoke specialized code. In all other cases we use the approach
% described below.
%
% For instance, we can use exponential distributions
%
% gnb('conditional','exponential') 
%
% or truncated gaussians
%
% truncgauss = @(x,mu,sigma)(normpdf(x,mu,sigma)./(normcdf(500,mu,sigma)-normcdf(0,mu,sigma)));
% truncmle = @(x)(mle(x,'pdf',truncgauss,'start', [nanmean(x) nanstd(x)],'lower', [0 0]));
% 
% gnb('conditional',truncgauss,'mle',truncmle)
% 
% In the former case, we can specify the conditional as a string (see code for more options). 
% In the latter case, both the conditional and maximum likelihood estimates are functions that 
% need to be specified. Check Matlab's mle function documentation for
% options.
%
% PARAMETERS:
%   params 
%   nclasses
%   priors
%
%   Copyright (c) 2008, Marcel van Gerven


    properties
      
      
      conditional = 'normal'; % existing pdf or function handle
      mle; % can be used to specify custom mle function
    
    end
    
    methods
      
      function obj = gnb(varargin)
        
        obj = obj@classifier(varargin{:});
      end
      
      function p = estimate(obj,X,Y)
        
        p.nclasses = obj.nunique(Y);
        nfeatures = size(X,2);
        
        % estimate class priors
        p.priors = zeros(p.nclasses,1);
        for j=1:p.nclasses
          p.priors(j) = sum(Y(:,1)==j)/size(Y,1);
        end
        
        % parameters per class/feature pair
        p.params = cell(p.nclasses,nfeatures);
        
        if ~isa(obj.conditional,'function_handle')
          
          % estimate class-conditional means
          for j=1:nfeatures
            for k=1:p.nclasses
              p.params{k,j} = mle(X(Y(:,1) == k,j),'distribution',lower(obj.conditional));
            end
          end
          
        else
          
          if isempty(obj.mle)
            
            for j=1:nfeatures
              for k=1:p.nclasses
                
                p.params{k,j} = mle(X(Y(:,1) == k,j),'pdf',obj.conditional);
              end
            end
            
          else
            
            for j=1:nfeatures
              for k=1:p.nclasses
                
                p.params{k,j} = obj.mle(X(Y(:,1) == k,j));
              end
            end
            
          end
        end
        
      end
      
      function Y = map(obj,X)
        
        p = obj.params.params;
        
        Y = nan(size(X,1),obj.params.nclasses);
        
        nparams = length(p{1,1});
        
        for m=1:size(Y,1) % iterate over examples
          
          for c=1:obj.params.nclasses
            
            % compute conditional
            conditional = zeros(1,size(X,2));
            
            if ~isa(obj.conditional,'function_handle')
              
              for j=1:size(X,2)
                if nparams == 1
                  conditional(j) = pdf(obj.conditional,X(m,j),p{c,j}(1));
                elseif nparams == 2
                  conditional(j) = pdf(obj.conditional,X(m,j),p{c,j}(1),p{c,j}(2));
                else
                  conditional(j) = pdf(obj.conditional,X(m,j),p{c,j}(1),p{c,j}(2),p{c,j}(3));
                end
              end
              
            else
              
              for j=1:size(X,2)
                if nparams == 1
                  conditional(j) = obj.conditional(X(m,j),p{c,j}(1));
                elseif nparams == 2
                  conditional(j) = obj.conditional(X(m,j),p{c,j}(1),p{c,j}(2));
                else
                  conditional(j) = obj.conditional(X(m,j),p{c,j}(1),p{c,j}(2),p{c,j}(3));
                end
              end
              
            end
            
            % degenerate cases
            if ~obj.params.priors(c) || any(isinf(conditional)) || ~all(conditional)
              Y(m,c) = 0;
              break
            end
            
            % compute probability
            Y(m,c) = log(obj.params.priors(c)) + mynansum(log(conditional));
            
          end
          
          % compute normalizing term using log-sum-exp trick
          
          mx = max(Y(m,:));
          
          nt = 0;
          for c=1:obj.params.nclasses
            nt = nt + exp(Y(m,c) - mx);
          end
          nt = log(nt) + mx;
          
          % normalize
          Y(m,:) = exp(Y(m,:) - nt);
          
        end
        
      end
      
    end
end
