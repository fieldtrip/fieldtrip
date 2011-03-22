classdef ft_mv_hmm < ft_mv_timeseries
%FT_MV_HMM Hidden Markov Model with discrete state and continuous observations
% 
% refs
% Pattern Recognition and Machine Learning, Bishop
% BNT toolbox
% 
% NOTE:
% - X and Y are represented as cell-arrays accommodates for multiple trials
% - X(t,k) denotes the t-th time point for the k-th observation
% - nhidden specifies the number of hidden states; this is used only when Y
%   is absent or NaN
%
% EXAMPLE:
% rand('seed',3); randn('seed',3);
% 
% nsamples = 1000; ncov = 10; ncycles = 10;
% Y = repmat([ones(100,1); 2*ones(100,1)],[5 1]);
% 
% k = ft_mv_hmm('verbose',true);
% X = repmat(Y,[1 ncov]); X(X(:)==1) = randn(1,sum(X(:)==1)); X(X(:)==2) = 0.1+3*randn(1,sum(X(:)==2));
% k = k.train(zscore(X),Y); % everything assumed observed
% X = repmat(Y,[1 ncov]); X(X(:)==1) = randn(1,sum(X(:)==1)); X(X(:)==2) = 0.1+3*randn(1,sum(X(:)==2));
% Z = k.predict(zscore(X));
% plot(Y,'k-');
% hold on;
% plot(Z,'ro');
% disp(mean(abs(Z - Y)));
%
% Copyright (c) 2010, Marcel van Gerven
  
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
    
    nhidden   % number of hidden states
    
    nmixture  % number of gaussian mixture components
    
  end

  methods
    
    function obj = ft_mv_hmm(varargin)
      
      obj = obj@ft_mv_timeseries(varargin{:});
      
      
    end

    function obj = train(obj,X,Y)

      if ~iscell(X), X = {X}; end
      for c=1:length(X)
        X{c} = X{c}';
      end
      nobs = size(X{1},1);
      
      if nargin > 2 
        
        if ~iscell(Y), Y = {Y}; end

        K = 0;
        for c=1:length(X)
          Y{c} = Y{c}';
          K = max(K,max(Y{c}));
        end

        obj.nhidden = 0;
        
      else        
        Y = {};
      end
      
      % in case of equal lengths we re-represent the timeseries
      % as features x time x repetitions; this is more efficient
      if length(unique(cellfun(@(x)(size(x,2)),X)))==1
        Xt = zeros([size(X{1}) length(X)]);
        for j=1:length(X)
          Xt(:,:,j) = X{j};
        end
        X = Xt;
        Y = cell2mat(Y);
      end

        
      if isempty(obj.nhidden) % Y assumed observed
        
        % class labels converted to state over time
        if size(Y,2) == 1
         Y = repmat(Y,[1 size(X,2)]);
        end

        nclasses = max(Y(:));
        
        covprior = repmat(obj.cov_prior*eye(nobs,nobs), [1 1 nclasses]);
        
        [obj.prior, obj.transmat, obj.mu, obj.Sigma] = gausshmm_train_observed(X, Y, nclasses,'cov_type',obj.cov_type,...
          'tied_cov',obj.tied_cov,'cov_prior',covprior);
        
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
          mhmm_em(X, prior0, transmat0, mu0, Sigma0, mixture,'cov_type',obj.cov_type,'max_iter',10);
        
      end
        
    end
        
    function post = test(obj,X)   
      % Viterbi decoding; outputs most likely path per sequence
     
      if ~iscell(X), X={X}; end
      
      post = cell(size(X));
      
      for c=1:size(X)
        
        if ~isempty(X{c})
          B = mixgauss_prob(X{c}', obj.mu, obj.Sigma);
        else
          B = ones(numel(obj.prior),size(X,2))/numel(obj.prior);
        end
        
        % NOTE: viterbi normalise clashes with afni normalise
        post{c} = viterbi_path(obj.prior, obj.transmat, B);
        
      end
      
    end
    
    function ll = likelihood(obj,X)
      % return the log likelihood of an observed sequence
      
      if ~iscell(X), X={X}; end
      
      ll = nan(length(X),1);
      
      for c=1:length(X)
      
         ll(c) = mhmm_logprob(X{c}', obj.prior, obj.transmat, obj.mu, obj.Sigma, obj.mixmat);

      end
      
    end
  
  end
  
end
