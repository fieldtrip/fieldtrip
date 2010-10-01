classdef ft_mv_naive < ft_mv_predictor
%FT_MV_NAIVE naive Bayes classifier; estimates discrete priors and
%class-conditional gaussian distributions for the individual parameters
%
% Reference:
% UPDATING FORMULAE AND A PAIRWISE ALGORITHM FOR COMPUTING SAMPLE VARIANCES
% by Tony F. Chan Gene H. Golub Randall J. LeVeque
%
% Note: sometimes we get complex values for sigma; to-be-solved; caused by
% instability of the online algorithm
% 
% Copyright (c) 2010, Marcel van Gerven

  properties
    
    n     % number of samples per class and feature
    S     % sum per class
    SS    % sum of squares per class
    
  end
  
  methods
    
    function obj = ft_mv_naive(varargin)

      obj = obj@ft_mv_predictor(varargin{:});

    end
    
    function obj = train(obj,X,Y)
  
      % multiple datasets
      if iscell(X) || iscell(Y)
        obj = ft_mv_ndata('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end  
     
      % multiple outputs
      if size(Y,2) > 1
        obj = ft_mv_noutput('mvmethod',obj);
        obj = obj.train(X,Y);
        return;
      end
        
      % at least two classes
      nclasses = max(2,max(Y));
      nfeatures = size(X,2);
   
      % initialize if not yet done
      if isempty(obj.n)
        obj.n = zeros(nclasses,nfeatures); 
        obj.S = zeros(nclasses,nfeatures);
        obj.SS = zeros(nclasses,nfeatures);
        binit = true;
      else
        if obj.verbose, fprintf('online learning\n'); end
        binit = false;
      end
      
      for k=1:nclasses
        
        % estimate class priors while taking nan into account
        obj.n(k,:) = obj.n(k,:) + sum(~isnan(X(Y==k,:)),1);
        
        % estimate class-conditional sum
        obj.S(k,:) = obj.S(k,:) + nansum(X(Y == k,:),1);
        
        % estimate class-conditional sum of squares
        mu    = obj.S ./ obj.n;
       
        if binit
        
          obj.SS(k,:) = sum(bsxfun(@minus,X(Y==k,:),mean(X(Y==k,:))).^2,1);
          
        else
          % use Eq 1.4 in UPDATING FORMULAE AND A PAIRWISE ALGORITHM FOR COMPUTING SAMPLE VARIANCES
          % by Tony F. Chan Gene H. Golub Randall J. LeVeque
          obj.SS(k,:) = obj.SS(k,:) + nansum(bsxfun(@minus,X(Y==k,:),mu(k,:)).^2,1) - ...
            bsxfun(@rdivide,(nansum(bsxfun(@minus,X(Y==k,:),mu(k,:)),1).^2),obj.n(k,:));
        end
        
      end
      
    end
    
    function Y = test(obj,X)

      prior = sum(obj.n,2)./sum(obj.n(:));

      % handle degenerate cases
      nzidx = find(prior ~= 0);

      % compute means and variances
      mu    = obj.S ./ obj.n;
      %sigma = (obj.SS ./ obj.n); % biased estimator
      sigma = obj.SS ./ (obj.n - 1);

      nclasses = length(prior);
      
      Y = zeros(size(X,1),nclasses);

      for m=1:size(Y,1) % iterate over examples

        for cc=1:length(nzidx)
          
          c = nzidx(cc);

          conditional = 1./sqrt(2*pi*sigma(c,:)) .* exp(- (X(m,:) - mu(c,:)).^2./(2*sigma(c,:)));

          % compute probability
          Y(m,c) = log(prior(c)) + nansum(log(conditional));

        end
        
        % compute normalizing term using log-sum-exp trick

        mx = max(Y(m,nzidx));

        nt = log(sum(exp(Y(m,nzidx)-mx))) + mx;

        % normalize
        Y(m,nzidx) = exp(Y(m,nzidx) - nt);

      end
      
    end

    function [m,desc] = model(obj)
      % return the parameters wrt a class label in some shape

      mu    = obj.S ./ obj.n;
      %sigma = (obj.SS ./ obj.n); % biased estimator
      sigma = obj.SS ./ (obj.n - 1);

      % return model for all classes; i.e., their means and standard
      % deviations
      if size(obj.n,1) == 2
        
        m(3:4) = mat2cell(mu,ones(1,size(mu,1)),size(mu,2));
        m(5:6) = mat2cell(sigma,ones(1,size(sigma,1)),size(sigma,2));
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
        desc{5} = 'variance for condition 1';
        desc{6} = 'variance for condition 2';
        
      else
        m = {};
        desc = {};
      end
        
    end

  end
  
end
