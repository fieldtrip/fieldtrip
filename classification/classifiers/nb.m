classdef nb < classifier
  %NB gaussian naive Bayes classifier; estimates discrete priors and
  %class-conditional gaussian distributions for the individual parameters
  %
  %   Copyright (c) 2008, Marcel van Gerven
  %
  %   $Log: nb.m,v $
  %

  properties

    priors
    means
    stds

  end

  methods
    function obj = nb(varargin)

      obj = obj@classifier(varargin{:});

    end
    function obj = train(obj,data,design)

      [data,design] = obj.check_input(data,design);
      
      if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end

      nfeatures = size(data,2);

      % estimate class priors
      obj.priors = zeros(obj.nclasses,1);
      for j=1:obj.nclasses
        obj.priors(j) = sum(design(:,1)==j)/size(design,1);
      end

      % estimate class-conditional means
      obj.means = zeros(obj.nclasses,nfeatures);
      for j=1:nfeatures
        for k=1:obj.nclasses
          obj.means(k,j) = nanmean(data(design(:,1) == k,j));
        end
      end

      % estimate class-conditional standard deviation
      obj.stds = zeros(obj.nclasses,nfeatures);
      for j=1:nfeatures
        for k=1:obj.nclasses
          obj.stds(k,j) = mynanstd(data(design(:,1) == k,j));
        end
      end

      % BUGFIX: handle zero SD
      obj.stds(obj.stds == 0) = 1e-20;

    end
    function post = test(obj,data)

      data = obj.check_input(data);
      
      post = nan(size(data,1),obj.nclasses);

      for m=1:size(post,1) % iterate over examples

        for c=1:obj.nclasses

          conditional = 1./(sqrt(2*pi)*obj.stds(c,:)) .* exp(- (data(m,:) - obj.means(c,:)).^2./(2*obj.stds(c,:).^2));

          % degenerate cases
          if ~obj.priors(c) || any(isinf(conditional)) || ~all(conditional)
            post(m,c) = 0;
            break
          end

          % compute probability
          post(m,c) = log(obj.priors(c)) + nansum(log(conditional));

        end

        % compute normalizing term using log-sum-exp trick

        mx = max(post(m,:));

        nt = 0;
        for c=1:obj.nclasses
          nt = nt + exp(post(m,c) - mx);
        end
        nt = log(nt) + mx;

        % normalize
        post(m,:) = exp(post(m,:) - nt);

      end
    end

    function m = getmodel(obj,label,dims)
      % return the parameters wrt a class label in some shape

      if nargin < 2 || isempty(label) || isnan(label)

        % return model for all classes; i.e., their means
        m = obj.means;

      else

        m = obj.means(label,:);

      end

      if numel(m) == prod(dims)
        m = reshape(m,dims);
      end

    end

  end
  
end
