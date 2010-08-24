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

      % return model for all classes; i.e., their means and standard
      % deviations
      if obj.params.nclasses == 2
        
        m(3:4) = mat2cell(obj.params.means,ones(1,size(obj.params.means,1)),size(obj.params.means,2));
        m(5:6) = mat2cell(obj.params.stds,ones(1,size(obj.params.stds,1)),size(obj.params.stds,2));
        m{1}   = 2*(m{3} - m{4}).^2 ./  (m{5}.^2 + m{6}.^2) + ...
          (m{5}.^2./m{6}.^2 - 1 - log(m{5}.^2./m{6}.^2))/2 + ...
          (m{6}.^2./m{5}.^2 - 1 - log(m{6}.^2./m{5}.^2))/2;
        
        C1 = m{5};
        C2 = m{6};
        C=(C1+C2)/2;
        dmu=(m{3}-m{4});
        m{2}=0.125*dmu.*C.*dmu + 0.5*log(C ./ sqrt(C1.*C2));
        
        desc = cell(1,5);
        desc{1} = 'Bhattacharyya distance'; % distance measure
        desc{2} = 'symmetricized KL divergence'; % Bhattacharyya distance alternative
        desc{3} = 'means for condition 1';
        desc{4} = 'means for condition 2';
        desc{5} = 'SD for condition 1';
        desc{6} = 'SD for condition 2';
        
      else
        m = {};
        desc = {};
      end
        
    end

  end
  
end
