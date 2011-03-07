classdef ft_mv_hmm < ft_mv_timeseries
%FT_MV_HMM Hidden Markov Model with discrete state and continuous observations
% 
% refs
% Pattern Recognition and Machine Learning, Bishop
% BNT toolbox
% 
% NOTE:
% X and Y as cell-arrays accommodates for multiple trials
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
     
    hmm; % learned hmm(s)
    
%     prior     % class prior
%     transmat  % transition matrix
%     mu        % conditional means
%     Sigma     % conditional covariance
    
    % if hidden > 0 then we learn a separate hmm for each class with nhidden
    % states. If nhidden = 0 then the state is assumed to be given by the
    % class variable such that the HMM is fully observed during training
    nhidden = 0;
    
  end

  methods
    
    function obj = ft_mv_hmm(varargin)
      
      obj = obj@ft_mv_timeseries(varargin{:});
      
      
    end

    function obj = train(obj,X,Y)

      % flip dims
      if iscell(X)
        K = 0;
        for c=1:length(X)
          sz = size(X);
          X{c} = X{c}';
          Y{c} = Y{c}';
          K = max(K,max(Y{c}));
        end
        nobs = size(X{1},1);
        
        % in case of equal lengths we re-represent the timeseries
        % as features x time x repetitions
        if length(unique(cellfun(@(x)(size(x,2)),X)))==1
          Xt = zeros([size(X{1}) length(X)]);
          for j=1:length(X)
            Xt(:,:,j) = X{j};
          end
          X = Xt;
          Y = cell2mat(Y);
        end
        
      else
        sz = size(X);
        X = reshape(X,[sz(1) prod(sz(2:end))])';
        Y = Y';
        nobs = size(X,1);
        K = max(Y);
      end
        
      if obj.nhidden == 0
        
        % class labels converted to state over time
        if size(Y,2) == 1
         Y = repmat(Y,[1 size(X,2)]);
        end
        
        covprior = repmat(obj.cov_prior*eye(nobs,nobs), [1 1 nclasses]);
        
        [obj.hmm.prior, obj.hmm.transmat, obj.hmm.mu, obj.hmm.Sigma] = gausshmm_train_observed(X, Y, nclasses,'cov_type',obj.cov_type,...
          'tied_cov',obj.tied_cov,'cov_prior',covprior);
        
      else

        % this could also be implemented with mixtures of gaussians
        % see: http://www.cs.ubc.ca/~murphyk/Software/HMM/hmm_usage.html
        
        Q = obj.nhidden; % number of hidden states
        
        obj.hmm = cell(1,K);
        for j=1:K
          
          X1 = X(:,:,Y==j);
          
          prior0 = normalise(rand(Q,1));
          transmat0 = mk_stochastic(rand(Q,Q));
          
          sz = size(X1); tmp = reshape(X1,[sz(1) prod(sz(2:end))]);
          mu0 = repmat(mean(tmp,2), [1 Q]);
          Sigma0 = repmat(cov(tmp'), [1 1 Q]);
          
          [obj.hmm{j}.LL, obj.hmm{j}.prior, obj.hmm{j}.transmat, obj.hmm{j}.mu, obj.hmm{j}.Sigma obj.hmm{j}.mixmat] = ...
            mhmm_em(X1, prior0, transmat0, mu0, Sigma0, [],'cov_type',obj.cov_type,'tied_cov');
        end
        
        
      end
        
    end
        
    function post = test(obj,X)   
      % Viterbi decoding
      
      % flip dims
      if iscell(X)
        post = cell(size(X));
        for c=1:length(X)
          post{c} = obj.test(X{c});
        end
        return
      else
        sz = size(X);
        X = reshape(X,[sz(1) prod(sz(2:end))])';
      end
        
      if obj.nhidden == 0
        
        if ~isempty(X)
          B = mixgauss_prob(X, obj.hmm.mu, obj.hmm.Sigma);
        else
          B = ones(numel(obj.hmm.prior),size(X,2))/numel(obj.hmm.prior);
        end
        
        % NOTE: viterbi normalise clashes with afni normalise
        res = viterbi_path(obj.hmm.prior, obj.hmm.transmat, B);
        res = res';
        
        K = numel(obj.hmm.prior); % number of states
        
        post = zeros(size(res,1),K);
        
        for c=1:K
          post(res == c,c) = 1;
        end
        
      else
        
        post = zeros(1,length(obj.hmm));
        for j=1:length(obj.hmm)
          post(j) = mhmm_logprob(X, obj.hmm{j}.prior, obj.hmm{j}.transmat, obj.hmm{j}.mu, obj.hmm{j}.Sigma, obj.hmm{j}.mixmat);
        end      
        post = (post == max(post));
        
      end
      
    end
  
  end
  
end
