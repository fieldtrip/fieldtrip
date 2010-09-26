classdef ft_mv_hmm < ft_mv_timeseries
%FT_MV_HMM Hidden Markov Model with discrete state and continuous observations
% 
% refs
% Pattern Recognition and Machine Learning, Bishop
% BNT toolbox
% 
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

  end

  methods
    
    function obj = ft_mv_hmm(varargin)
      
      obj = obj@ft_mv_timeseries(varargin{:});
      
      
    end

    function obj = train(obj,X,Y)
              
      % flip
      if iscell(X)
        nclasses = 0;
        for c=1:length(X)
          X{c} = X{c}';
          Y{c} = Y{c}';
          nclasses = max(nclasses,max(Y{1}));        
        end
        nobs = size(X{1},1);        
      else
        X = X';
        Y = Y';
        nobs = size(X,1);
        nclasses = max(Y);
      end
      
      cov_prior = repmat(obj.cov_prior*eye(nobs,nobs), [1 1 nclasses]);
      [obj.prior, obj.transmat, obj.mu, obj.Sigma] = gausshmm_train_observed(X, Y, nclasses,'cov_type',obj.cov_type,...
        'tied_cov',obj.tied_cov,'cov_prior',cov_prior);
            
    end
        
    function post = test(obj,X)   
      % Viterbi decoding
      
      % flip
      if iscell(X)
        for c=1:length(X)
          X{c} = X{c}';
        end
      else
        X = X';
      end
      
      if ~isempty(X)      
        B = mixgauss_prob(X, obj.mu, obj.Sigma);
      else
        B = ones(numel(obj.prior),size(X,2))/numel(obj.prior);
      end
      res = viterbi_path(obj.prior, obj.transmat, B);
      res = res';
      
      post = zeros(size(res,1),2);
      
      for c=1:2
        post(res == c,c) = 1;
      end
      
    end
  
  
  end
  
end
