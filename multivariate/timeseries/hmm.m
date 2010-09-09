classdef hmm < classifier & timeseries
%HMM Hidden Markov Model with discrete state and continuous observations
% 
%
% refs
% Pattern Recognition and Machine Learning, Bishop
% BNT toolbox
%
%
% Copyright (c) 2010, Marcel van Gerven
  
  properties
    
    prior; % used to override learnt initial prior on states
    
    cov_type = 'diag'; % covariance type (full, diag, spherical)
    tied_cov = 1; % tie the covariances

    nclasses = 2; % number of hidden states
    
  end

  methods
    
    function obj = hmm(varargin)
      
      obj = obj@classifier(varargin{:});
      
      
    end

    function p = estimate(obj,X,Y)
              
      % flip
      if iscell(X)
        for c=1:length(X)
          X{c} = X{c}';
          Y{c} = Y{c}';
        end
      else
        X = X';
        Y = Y';
      end
      
      [p.prior, p.transmat, p.mu, p.Sigma] = gausshmm_train_observed(X, Y, 2,'cov_type',obj.cov_type,'tied_cov',obj.tied_cov);
      
      if ~isempty(obj.prior)
        p.prior = obj.prior;
      end
      
    end
        
    function post = map(obj,X)   
      % Viterbi decoding
      
      % flip
      if iscell(X)
        for c=1:length(X)
          X{c} = X{c}';
        end
      else
        X = X';
      end
      
      B = mixgauss_prob(X, obj.params.mu, obj.params.Sigma);
      res = viterbi_path(obj.params.prior, obj.params.transmat, B);
      
      res = res';
      
      post = zeros(size(res,1),2);
      
      for c=1:2
        post(res == c,c) = 1;
      end
      
    end
  
  
  end
  
end
