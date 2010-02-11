classdef hgnb_transfer < classifier & transfer_learner
  %hierarchical gaussian naive Bayes classifier
  %
  %   This classifier can be used for transfer learning:
  %
  %   y_si    ~ N(theta_s,sigma^2)
  %   theta_s ~ N(mu, tau^2)
  %
  %   Gelman, Bayesian data analysis, pp. 134
  %
  %   Copyright (c) 2009, Marcel van Gerven
  

  properties

    priors
    examples
    means
    stds
    
    nclasses

  end

  methods
    
    function obj = hgnb_transfer(varargin)

      obj = obj@classifier(varargin{:});

    end
    
    function p = estimate(obj,X,Y)

      p.nclasses = obj.nunique(Y{1});
      
      ntasks = length(X);

      nfeatures = size(X{1},2);

      p.means = cell(1,ntasks);
      p.stds  = cell(1,ntasks);

      for c=1:ntasks
       
        p.means{c} = zeros(p.nclasses,nfeatures);
        p.stds{c} = zeros(p.nclasses,nfeatures);
      end

      % estimate class priors by pooling over subjects
      p.priors = zeros(p.nclasses,1);
      for k=1:p.nclasses
        for c=1:ntasks
          p.priors(k) = p.priors(k) + sum(Y{c}(:,1)==k);
        end
      end
      p.priors = p.priors ./ sum(p.priors);

      alldata   = cat(1,X{:});
      alldesign = cat(1,Y{:});
      
      % iterate over features and classes
      for j=1:nfeatures

        for k=1:p.nclasses

          % compute per subject examples
          ns = zeros(1,ntasks);
          for c=1:ntasks
            ns(c) = length(Y{c}(:,1));
          end
          n = sum(ns);
          
          % compute variance by pooling over subjects
          sigma2 = mynanstd(alldata(alldesign(:,1) == k,j)).^2;
                
          % compute per subject variances
          sigmas2 = sigma2./ns;
%           sigmas2 = zeros(1,ntasks);
%           for c=1:ntasks
%             sigmas2(c) = mynanstd(data{c}(design{c}(:,1) == k,j)).^2;
%           end

          % compute per subject sample means
          ys = zeros(1,ntasks);
          for c=1:ntasks
            ys(c) = mynanmean(X{c}(Y{c}(:,1) == k,j));
          end

          % compute pooled estimate y..
          % acts as point estimate for hyperparameter mu
          y = sum(ys./sigmas2) / sum(1./sigmas2);

          % point estimate of hyperparameter tau using MSB and MSW

          MSB = sum((ys - y).^2)/(ntasks-1);

          MSW = 0;
          for c=1:ntasks
            MSW = MSW + sum((X{c}(Y{c}(1,:) == k,j) - ys(c)).^2);
          end
          MSW = MSW / (ntasks*(n-1));

          tau2 = MSB - MSW / n;

          for c=1:ntasks

            p.means{c}(k,j) = (ys(c)./sigmas2(c) + y/tau2) / (1./sigmas2(c) + 1/tau2);

%             obj.stds{c}(k,j) = sqrt(sigmas2(c));
            p.stds{c}(k,j) = sqrt(sigma2);

            % BUGFIX: handle zero SD
            p.stds{c}(p.stds{c} ==0) = 1e-20;
          end

        end
      end
    end

    function Y = map(obj,X)

      ntasks = length(X);

      Y = cell(1,ntasks);
      for c=1:ntasks

        Y{c} = nan(size(X{c},1),obj.params.nclasses);

        for m=1:size(Y{c},1) % iterate over examples

          for k=1:obj.params.nclasses

            conditional = 1./(sqrt(2*pi)*obj.params.stds{c}(k,:)) .* exp(- (X{c}(m,:) ...
              - obj.params.means{c}(k,:)).^2./(2*obj.params.stds{c}(k,:).^2));

            % degenerate cases
            if ~obj.params.priors(k) || any(isinf(conditional)) || ~all(conditional)
              Y{c}(m,k) = 0;
              break
            end

            % compute probability
            Y{c}(m,k) = log(obj.params.priors(k)) + mynansum(log(conditional));

          end

          % compute normalizing term using log-sum-exp trick

          mx = max(Y{c}(m,:));

          nt = 0;
          for k=1:obj.nclasses
            nt = nt + exp(Y{c}(m,k) - mx);
          end
          nt = log(nt) + mx;

          % normalize
          Y{c}(m,:) = exp(Y{c}(m,:) - nt);

        end
        
        
      end

    end
    
    function [m,desc] = getmodel(obj)
      % return the parameters
        
        m = cell(size(obj.params.means{1},1),length(obj.params.means));
        for c=1:size(m,1)
          for j=1:size(m,2)
            m{c,j} = full(obj.params.means{j}(c,:)); % ignore bias term
          end
        end
        
        desc = {'unknown'};
           
    end

  end
end
