classdef nb < classifier
  %NB gaussian naive Bayes classifier; estimates discrete priors and
  %class-conditional gaussian distributions for the individual parameters
  %
  %   Copyright (c) 2008, Marcel van Gerven

  methods
    
    function obj = nb(varargin)

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

      % estimate class-conditional means
      p.means = zeros(p.nclasses,nfeatures);
      for j=1:nfeatures
        for k=1:p.nclasses
          p.means(k,j) = mynanmean(X(Y(:,1) == k,j));
        end
      end

      % estimate class-conditional standard deviation
      p.stds = zeros(p.nclasses,nfeatures);
      for j=1:nfeatures
        for k=1:p.nclasses
          p.stds(k,j) = mynanstd(X(Y(:,1) == k,j));
        end
      end

      % BUGFIX: handle zero SD
      p.stds(p.stds == 0) = 1e-20;

    end
    
    function Y = map(obj,X)

       Y = nan(size(X,1),obj.params.nclasses);

      for m=1:size(Y,1) % iterate over examples

        for c=1:obj.params.nclasses

          conditional = 1./(sqrt(2*pi)*obj.params.stds(c,:)) .* exp(- (X(m,:) - obj.params.means(c,:)).^2./(2*obj.params.stds(c,:).^2));

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

    function [m,desc] = getmodel(obj)
      % return the parameters wrt a class label in some shape

      % return model for all classes; i.e., their means
        m = mat2cell(obj.params.means,ones(1,size(obj.params.means,1)),size(obj.params.means,2));
        desc = repmat({'unknown'},size(m));
        
    end

  end
  
end
