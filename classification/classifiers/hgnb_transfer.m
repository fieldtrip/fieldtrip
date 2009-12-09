classdef hgnb_transfer < transfer_classifier
  %HGNB hierarchical gaussian naive Bayes classifier
  %
  %   This classifier can be used for transfer learning:
  %
  %   y_si    ~ N(theta_s,sigma^2)
  %   theta_s ~ N(mu, tau^2)
  %
  %   Gelman, Bayesian data analysis, pp. 134
  %
  %   Copyright (c) 2009, Marcel van Gerven
  %
  %   $Log: hgnb.m,v $
  %

  properties

    priors
    examples
    means
    stds

  end

  methods
    function obj = hgnb_transfer(varargin)

      obj = obj@classifier(varargin{:});

    end
    function obj = train(obj,data,design)

      [data,design] = obj.check_input(data,design);
      
      ntasks = length(data);

      if isnan(obj.nclasses), obj.nclasses = max(design{1}(:,1)); end

      nfeatures = size(data{1},2);

      obj.means = cell(1,ntasks);
      obj.stds  = cell(1,ntasks);

      for c=1:ntasks
        obj.means{c} = zeros(obj.nclasses,nfeatures);
        obj.stds{c} = zeros(obj.nclasses,nfeatures);
      end

      % estimate class priors by pooling over subjects
      obj.priors = zeros(obj.nclasses,1);
      for k=1:obj.nclasses
        for c=1:ntasks
          obj.priors(k) = obj.priors(k) + sum(design{c}(:,1)==k);
        end
      end
      obj.priors = obj.priors ./ sum(obj.priors);

      alldata   = cat(1,data{:});
      alldesign = cat(1,design{:});
      
      % iterate over features and classes
      for j=1:nfeatures

        for k=1:obj.nclasses

          % compute per subject examples
          ns = zeros(1,ntasks);
          for c=1:ntasks
            ns(c) = length(design{c}(:,1));
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
            ys(c) = mynanmean(data{c}(design{c}(:,1) == k,j));
          end

          % compute pooled estimate y..
          % acts as point estimate for hyperparameter mu
          y = sum(ys./sigmas2) / sum(1./sigmas2);

          % point estimate of hyperparameter tau using MSB and MSW

          MSB = sum((ys - y).^2)/(ntasks-1);

          MSW = 0;
          for c=1:ntasks
            MSW = MSW + sum((data{c}(design{c}(1,:) == k,j) - ys(c)).^2);
          end
          MSW = MSW / (ntasks*(n-1));

          tau2 = MSB - MSW / n;

          for c=1:ntasks

            obj.means{c}(k,j) = (ys(c)./sigmas2(c) + y/tau2) / (1./sigmas2(c) + 1/tau2);

%             obj.stds{c}(k,j) = sqrt(sigmas2(c));
            obj.stds{c}(k,j) = sqrt(sigma2);

            % BUGFIX: handle zero SD
            obj.stds{c}(obj.stds{c} ==0) = 1e-20;
          end

        end
      end
    end

    function post = test(obj,data)

      data = obj.check_input(data);
      
      ntasks = length(data);

      post = cell(1,ntasks);
      for c=1:ntasks

        post{c} = nan(size(data{c},1),obj.nclasses);

        for m=1:size(post{c},1) % iterate over examples

          for k=1:obj.nclasses

            conditional = 1./(sqrt(2*pi)*obj.stds{c}(k,:)) .* exp(- (data{c}(m,:) ...
              - obj.means{c}(k,:)).^2./(2*obj.stds{c}(k,:).^2));

            % degenerate cases
            if ~obj.priors(k) || any(isinf(conditional)) || ~all(conditional)
              post{c}(m,k) = 0;
              break
            end

            % compute probability
            post{c}(m,k) = log(obj.priors(k)) + mynansum(log(conditional));

          end

          % compute normalizing term using log-sum-exp trick

          mx = max(post{c}(m,:));

          nt = 0;
          for k=1:obj.nclasses
            nt = nt + exp(post{c}(m,k) - mx);
          end
          nt = log(nt) + mx;

          % normalize
          post{c}(m,:) = exp(post{c}(m,:) - nt);

        end
      end

    end
    
    function m = getmodel(obj,label,dims)
      % return the parameters wrt a class label in some shape

      if nargin < 2 || isempty(label) || isnan(label)

        % return model for all classes; i.e., their means
        m = obj.means;

      else

        m = obj.means;
        for c=1:length(m)
          m{c} = m{c}(label,:);
        end

      end      

      for c=1:length(m)
        if numel(m) == prod(dims)
          m{c} = reshape(m,dims);
        end
      end
      
    end

  end
end
