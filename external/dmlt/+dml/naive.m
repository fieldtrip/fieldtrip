classdef naive < dml.method
% NAIVE gaussian naive Bayes classifier.
%
%   DESCRIPTION
%   Estimates discrete priors and class-conditional gaussian distributions 
%   for the individual parameters.
%
%   Suitable for online learning (repeated calls to the train function).
%   Note: sometimes we get complex values for sigma; to-be-solved; caused by
%   instability of the online algorithm
%
%   REFERENCE
%   Updating formulae and a pairwise algorithm for computing sample
%   variances by Tony F. Chan Gene H. Golub Randall J. LeVeque
%
%   EXAMPLE
%   X = rand(10,20); Y = [1 1 1 1 1 2 2 2 2 2]';
%   m = dml.naive
%   m = m.train(X,Y);
%   Z = m.test(X);
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)


  properties
    
    n     % number of samples per class and feature
    S     % sum per class
    SS    % sum of squares per class
    
  end
  
  methods
    
    function obj = naive(varargin)

      obj = obj@dml.method(varargin{:});

    end
    
    function obj = train(obj,X,Y)
        
      % handle multiple datasets
      if iscell(X)
        obj = dml.ndata('method',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % multiple outputs
      if size(Y,2) > 1
        obj = dml.noutput('method',obj);
        obj = obj.train(X,Y);
        return;
      end
      
      % at least two classes
      nclasses = max(2,max(Y));
      nfeatures = size(X,2);
   
      % initialize if not yet done
      if obj.restart || isempty(obj.n)
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
        mu = obj.S(k,:) ./ obj.n(k,:);
       
        if binit
        
          obj.SS(k,:) = sum(bsxfun(@minus,X(Y==k,:),mean(X(Y==k,:))).^2,1);
          
        else
          % use Eq 1.4 in UPDATING FORMULAE AND A PAIRWISE ALGORITHM FOR COMPUTING SAMPLE VARIANCES
          % by Tony F. Chan Gene H. Golub Randall J. LeVeque
          obj.SS(k,:) = obj.SS(k,:) + nansum(bsxfun(@minus,X(Y==k,:),mu).^2,1) - ...
            bsxfun(@rdivide,(nansum(bsxfun(@minus,X(Y==k,:),mu),1).^2),obj.n(k,:));
        end
        
      end
      
    end
    
    function Y = test(obj,X)

      prior = sum(obj.n,2)./sum(obj.n(:));

      % handle degenerate cases
      nzidx = find(prior ~= 0);

      % compute means and variances
      mu    = obj.S ./ obj.n;
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

    function m = model(obj)
      % returns
      %
      % m.divergence summed KL divergences between all pairs of class distributions

      nclasses = size(obj.n,1);
      nfeatures = size(obj.n,2);
      
      mu    = obj.S ./ obj.n;
      sigma = obj.SS ./ (obj.n - 1);

      d = zeros(1,nfeatures);
      for i=1:nclasses
        for j=[1:(i-1) (i+1):nclasses]

          d = d + log(sigma(j,:)./sigma(i,:)) + (sigma(i,:).^2 + (mu(i,:)-mu(j,:)).^2) ./ (2*sigma(j,:).^2) - 0.5;
          
        end
      end
        
      m = [];
      m.divergence = d(:);
      
    end

  end
  
end
