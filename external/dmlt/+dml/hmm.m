classdef hmm < dml.method
% HMM Hidden Markov model with continuous observations.
% 
% DESCRIPTION
% input X is of size repetitions x features x timepoints
% note:
% - nhidden specifies the number of hidden states; this is used only when Y
%   is absent or NaN
%
% EXAMPLE
% rand('seed',3); randn('seed',3);
% 
% nsamples = 1000; ncov = 10; ncycles = 10;
% Y = repmat([ones(100,1); 2*ones(100,1)],[5 1]);
% 
% k = dml.hmm('verbose',true);
% X = repmat(Y,[1 ncov]); X(X(:)==1) = randn(1,sum(X(:)==1)); X(X(:)==2) = 0.1+3*randn(1,sum(X(:)==2));
% X = reshape(X',[1 size(X,2) size(X,1)]);
% Y = reshape(Y',[1 size(Y,2) size(Y,1)]);
% k = k.train(X,Y); % everything assumed observed
% Z = k.test(X);
% plot(squeeze(Y),'k-');
% hold on;
% plot(squeeze(Z),'ro');
%
% REFERENCES
% Pattern Recognition and Machine Learning, Bishop
% BNT toolbox
%
% DEVELOPER
% Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
    cov_type = 'full';     % covariance type (full, diag, spherical)
    tied_cov =  0;         % tie the covariances
    cov_prior = 0;         % prior on the covariance matrix (diagonal term)
    
    prior     % class prior
    transmat  % transition matrix
    mu        % conditional means
    Sigma     % conditional covariance
    mixmat    % mixing matrix
    LL        % likelihood
    
    nhidden = 2  % number of hidden states; state assumed observed if hidden
    
    nmixture  % number of gaussian mixture components
    
    maxiter = 100; % max number of em iterations
    
  end

  methods
    
    function obj = hmm(varargin)
      
      obj = obj@dml.method(varargin{:});
      
      
    end

    function obj = train(obj,X,Y)

      % data represented as repetitions x features x timepoints

      if ~isempty(obj.indims)
        X = reshape(X,[size(X,1) obj.indims]);
      end
            
      X = permute(X,[2 3 1]);
        
      if isempty(obj.nhidden) || (nargin==3 && (~all(isnan(Y(:)))))% Y assumed observed

        Y = permute(Y,[2 3 1]);

        nclasses = max(Y(:));
        nobs = size(X,1);
        
        covprior = repmat(obj.cov_prior*eye(nobs,nobs), [1 1 nclasses]);
        
        [obj.prior, obj.transmat, obj.mu, obj.Sigma] = gausshmm_train_observed(X, Y, nclasses,'cov_type',obj.cov_type,...
          'tied_cov',obj.tied_cov,'cov_prior',covprior);
        obj.mixmat = ones(size(obj.prior,1),1);
        
      else % Y assumed unobserved

        % this could also be implemented with mixtures of gaussians
        % see: http://www.cs.ubc.ca/~murphyk/Software/HMM/hmm_usage.html
        
        Q = obj.nhidden; % number of hidden states
                
        prior0 = normalise(rand(Q,1));
        transmat0 = mk_stochastic(rand(Q,Q));
        
        sz = size(X); tmp = reshape(X,[sz(1) prod(sz(2:end))]);
        if isempty(obj.nmixture)
          mu0 = repmat(mean(tmp,2), [1 Q]);
          Sigma0 = repmat(cov(tmp'), [1 1 Q]);
          mixture = [];
        else       
          % NOTE: mixture of gaussians requires NETLAB
          % http://www.mathworks.com/matlabcentral/fileexchange/2654-netlab
          O = size(X,1);
          [mu0, Sigma0] = mixgauss_init(Q*obj.nmixture, tmp, obj.cov_type);
          mu0 = reshape(mu0, [O Q obj.nmixture]);
          Sigma0 = reshape(Sigma0, [O O Q obj.nmixture]);
          mixture = mk_stochastic(rand(Q,obj.nmixture));
        end
        
        [obj.LL, obj.prior, obj.transmat, obj.mu, obj.Sigma obj.mixmat] = ...
          mhmm_em(X, prior0, transmat0, mu0, Sigma0, mixture,'cov_type',obj.cov_type,'max_iter',obj.maxiter);
        
      end
        
    end
        
    function post = test(obj,X)   
      % Viterbi decoding; outputs most likely path per sequence
     
      post =zeros(size(X,1),1,size(X,3));
      
      for c=1:size(X,1)
        
        B = mixgauss_prob(reshape(X(c,:,:),[size(X,2) size(X,3)]), obj.mu, obj.Sigma);
        
        % NOTE: viterbi normalise clashes with afni normalise
        post(c,:,:) = viterbi_path(obj.prior, obj.transmat, B);
        
      end
            
    end
    
    function ll = likelihood(obj,X)
      % return the log likelihood of an observed sequence
      
      ll = nan(size(X,1),1);
      
      for c=1:size(X,1)
      
         ll(c) = mhmm_logprob(reshape(X(c,:,:),[size(X,2) size(X,3)]), obj.prior, obj.transmat, obj.mu, obj.Sigma, obj.mixmat);

      end
      
    end
  
  end
  
end
