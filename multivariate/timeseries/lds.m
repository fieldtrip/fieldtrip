classdef lds < regressor & timeseries & transfer_learner
%LDS linear dynamical system 
%
% state can be fully observed/unobserved during training.
% partial observability of multiple observations is supported
% observations are normally distributed conditional on the state.
% multiple observations sequences are supported
%
% bias term is *NOT* automatically added to the model
%
%
% refs
% Pattern Recognition and Machine Learning, Bishop
% A unifying review of linear dynamical systems, Gharamani 
% IEEE: ...
%
% options:
%  maxiter = 100; % maximum number of EM iterations
%  thresh = 1e-4; % EM convergence threshold    
%  diagQ = 0; % regularize state noise to diagonal (0 <= diagQ <=1)
%  diagR = 0; % regularize measurement noise to diagonal (0 <= diagR <=1)
%
% parameters:
%
% K = number of states
% M = number of observations
% T = number of timesteps
%
% A     % K x K transition matrix for the unobservable state
% C     % M x K emission matrix for the observations
%
% mu0   % K x 1 the initial mean of the hidden state
% V0    % K x K the initial hidden noise covariance
% R     % M x M measurement noise covariance
% Q     % K x K state noise covariance
% 
% epsilon % added to the diagonal of the covariance matrices for numerical
%           stability
%
% NOTE: X and Y are swapped wrt the Kalman filter conventions
%
% EXAMPLE:
%
% rand('seed',1); randn('seed',1);
% nsamples = 1000; ncov = 10; ncycles = 10;
% Y = sin(ncycles * 2 * pi * (1:nsamples) ./ nsamples)';
% 
% k = lds('verbose',true);
% X = repmat(Y,[1 ncov]) + randn(size(Y,1),ncov); 
% k = k.train(zscore(X),zscore(Y)); % everything assumed observed
% X = repmat(Y,[1 ncov]) + randn(size(Y,1),ncov); 
% Z = k.test(zscore(X));
% plot(Z,'ko');
% disp(mean(abs(Z - Y)));
% 
% k = lds('verbose',true);
% k = k.train(zscore(X),nan(size(Y))); % hidden state estimation
% X = repmat(Y,[1 ncov]) + randn(size(Y,1),ncov); 
% plot(k.test(zscore(X)),'ko');
% 
%
% Copyright (c) 2010, Marcel van Gerven
  
  properties
    
    maxiter = 100; % maximum number of EM iterations
    thresh = 1e-4; % EM convergence threshold    
    diagQ = 0; % diagonal state noise
    diagR = 0; % diagonal measurement noise
    
    inference = 'smooth'; % filter / smooth
    
    epsilon = 1e-7;
    
  end

  methods
    
    function obj = lds(varargin)
      
      obj = obj@regressor(varargin{:});
    end
    
    function p = estimate(obj,X,Y)
               
      % cast to cell array (multiple observation sequences)
      % representation as variables * timepoints
      if ~iscell(X)
        X = {X'}; 
      else
        for c=1:length(X)
          X{c} = X{c}';
        end
      end
      
      if nargin < 3
        
        % assume one signal underlying the observations
        Y = cell(1,length(X));
        for c=1:length(X)
          Y{c} = nan(1,size(X,2));
        end
        
      else
        
        if ~iscell(Y)
          Y = {Y'};
        else
          for c=1:length(Y)
            Y{c} = Y{c}';
          end
        end
           
      end
      
      % remove empty elements
      eidx = cellfun(@(x)(isempty(x)),X);
      X = X(~eidx);
      Y = Y(~eidx);
        
      K = size(Y{1},1);
      M = size(X{1},1);
      N = length(X);
    
      if isempty(obj.params)
        % initialize to small random values
        
        obj.params.A    = 1e-3*randn(K,K);
        obj.params.C    = 1e-3*randn(M,K);
        obj.params.mu0  = 1e-3*randn(K,1);
        obj.params.V0   = 1e-3*rand(K,K);
        obj.params.R    = 1e-3*rand(M,M);
        obj.params.Q    = 1e-3*rand(K,K);
        
        % symmetricize 
        obj.params.V0 = (obj.params.V0 + obj.params.V0') ./ 2;
        obj.params.R = (obj.params.R + obj.params.R') ./ 2;
        obj.params.Q = (obj.params.Q + obj.params.Q') ./ 2;
        
      end
      
      oldLL = -inf;
      LL = 0;
      loglik = [];
      iter = 0;
      while abs(LL - oldLL) > obj.thresh && iter < obj.maxiter
        
        oldLL = LL;
        
        % E step
        LL=0;
        for c=1:N % iterate over sequences
          
          [G1t,G2t,G3t,G4t,G5t,G6t,mu0t,V0t,LLt] = obj.Estep(X{c},Y{c});
          
          if c==1
          
            G1=G1t; 
            G2=G2t; 
            G3=G3t; 
            G4=G4t; 
            G5=G5t; 
            G6=G6t; 
            mu0=mu0t; 
            V0=V0t + mu0t*mu0t'; 
            LL=LLt;
          
          else
            
            G1 = G1 + G1t;
            G2 = G2 + G2t;
            G3 = G3 + G3t;
            G4 = G4 + G4t;
            G5 = G5 + G5t;
            G6 = G6 + G6t;
            
            mu0 = mu0 + mu0t;
            V0 = V0 + V0t + mu0t*mu0t';

            LL = LL + LLt;
            
          end
               
        end
        
        loglik = [loglik LL];
        
        if obj.verbose && ~isnan(LL)
          fprintf('EM step: %d; loglikelihood: %g\n',iter,LL);
        end
        
        % M step
        
        T = sum(cellfun(@(x)(size(x,2)),X));
        
        p.C = G6 / G1;
        
        p.R = (G5 - p.C * G6') ./ T;
        p.R = (p.R + p.R') ./ 2;
        
        if obj.diagR
          p.R = (1-obj.diagR) * p.R + obj.diagR * diag(diag(p.R));
        end
        
        p.A = G4 / G3;
        
        p.Q = (G2 - p.A * G4') ./ (T-N);
        p.Q = (p.Q + p.Q') ./ 2;
        
        if obj.diagQ
          p.Q = (1-obj.diagQ) * p.Q + obj.diagQ * diag(diag(p.Q));
        end
        
        p.mu0 = mu0/N;
        p.V0 = V0/N - p.mu0*p.mu0';      
        
        % symmetricize and make positive semidefinite
        p.V0 = (p.V0 + p.V0') ./ 2;
        p.V0 = p.V0 + obj.epsilon*eye(size(p.V0));
     
        obj.params = p;
        
        iter = iter + 1;
        
        if ~any(isnan(Y{1}(:)))
          break;
        end
        
      end
    
      LL=0;
      for c=1:N
        [mu1,V1] = obj.filter(X{c});
        [mu,V,J,LLt] = obj.smooth(X{c},mu1,V1);
        LL = LL+LLt;
      end
      loglik = [loglik LL];
      
      if obj.verbose
        fprintf('EM step: %d; loglikelihood: %g\n',iter,LL);
      end
      
      p.loglik = loglik(2:end);
      
      
    end
    
    function [mu,V,loglik] = map(obj,X)   
     % LDS inference

     if ~iscell(X), X = {X}; end
     
     % remove empty elements
     neidx = cellfun(@(x)(~isempty(x)),X);
     idx = find(neidx);
     X = X(idx);
      
     mu = cell(1,length(neidx));
     for c=1:length(X)
     
       % representation as variables * timepoints
       XX = X{c}';
       
       [m,V] = obj.filter(XX);
       
       if strcmp(obj.inference,'smooth')
         [m,V,J,loglik] = obj.smooth(XX,m,V);
       end
       
       mu{idx(c)} = m';
     
     end
     
     if length(neidx)==1
       mu = mu{1};
     end
     
    end
    
    function [mu,V,J,loglik] = smooth(obj,X,mu1,V1)
      % Kalman smoother; uses filtered means and variances
      
      T = numel(V1);
      
      A = obj.params.A;
      V0 = obj.params.V0;
      Q = obj.params.Q;
      
      % RTS equations
      
      mu = mu1;
      V = V1;
      J = cell(1,T-1); % needed to compute transition matrix A in EM
      for n=(T-1):-1:1
        
        if n==1
          P = V0;
        else
          P = A * V1{n} * A' + Q;
        end
        
        if ~all(V{n}(:)==0)
          
          J{n} = (V1{n} * A') / P;
          
          mu(:,n) = mu1(:,n) + J{n} * (mu(:,n+1) - A * mu1(:,n));
          V{n} = V1{n} + J{n} * (V{n+1} - P)*J{n}';
          
        else
          
          % fix when variance is zero (deterministic state)
          J{n} = zeros(size(V{n}));
          mu(:,n) = mu1(:,n);
          V{n} = zeros(size(V{n}));
          
        end
        
        % symmetricize and make positive semidefinite
        V{n} = (V{n} + V{n}') ./ 2;
        V{n} = V{n} + obj.epsilon*eye(size(V{n}));

        
      end
      
      if nargout >= 4
        loglik = obj.compute_loglik(X,mu,V);
      end
      
    end
    
    function [mu,V,LL] = filter(obj,X) %HERE!!!!!! + first estimate params!   
      % Kalman filter algorithm; returns mean mu and covariance matrix V of
      % the states, the log likelihood and the observations Xobs with missing
      % values imputed (the latter remains to be done)
      
      A = obj.params.A;
      C = obj.params.C;
      mu0 = obj.params.mu0;
      V0 = obj.params.V0;
      R = obj.params.R;
      Q = obj.params.Q;

      K = size(A,1);
      T = size(X,2);
      
      mu = nan(K,T);
      V = cell(1,T);
      I = eye(K);      

      % take unobserved measurements into account

      obs = ~isnan(X(:,1));
      Cb = C(obs,:); % emission matrix for observed measurements
      
      % Kalman gain matrix
      K = (V0 * Cb') / (Cb * V0 * Cb' + R(obs,obs));

      mu(:,1) = mu0 + K * (X(obs,1) - Cb * mu0);
      V{1}  = (I - K * Cb) * V0;
      
      % symmetricize and make positive semidefinite
      V{1} = (V{1} + V{1}') ./ 2;
      V{1} = V{1} + obj.epsilon*eye(size(V{1}));
     
      for t=2:T

        obs = ~isnan(X(:,t));  
        Cb = C(obs,:); % emission matrix for observed measurements
        
        P = A * V{t-1} * A' + Q;
        
        PC = P * Cb';
      
        K = PC / (Cb * PC + R(obs,obs));
        
        AM = A * mu(:,t-1);
        mu(:,t) = AM + K * (X(obs,t) - Cb * AM);
        V{t}  = (I - K * Cb) * P;
     
        % symmetricize and make positive semidefinite
        V{t} = (V{t} + V{t}') ./ 2;
        V{t} = V{t} + obj.epsilon*eye(size(V{t}));
  
      end
      
      % compute log likelihood using error decomposition ignoring constant terms  
      % why don't we need the smoothed estimates of mu? CHECK THIS
      if nargout >= 3
        LL = obj.compute_loglik(X,mu,V);
      end
      
    end
        
  end
  
  methods(Access=protected)
    
    function [G1,G2,G3,G4,G5,G6,mu0,V0,LL] = Estep(obj,X,Y)
              
      K = size(Y,1);
      M = size(X,1);
      T = size(X,2);
      
      if any(isnan(Y(:))) % hidden state
                
        [mu1,V1] = obj.filter(X);
        [Y,V,J,LL] = obj.smooth(X,mu1,V1);

      else
      
        V = repmat({0},size(Y));
        J = repmat({0},[size(Y,1) size(Y,2)-1]);
        LL = nan;
        
      end
      
      G1 = zeros(K,K);
      G4 = zeros(K,K);
      for t=1:T
        
        G1  = G1 + V{t} + Y(:,t)*Y(:,t)';
        
        if t>1 % assuming Bishop's notation is correct
          G4 = G4 + J{t-1}*V{t} + Y(:,t)*Y(:,t-1)';
        end
        
      end
      
      G2 = G1 - Y(:,1) * Y(:,1)' - V{1};
      G3 = G1 - Y(:,T) * Y(:,T)' - V{T};
      
      mu0 = Y(:,1);
      V0 = V{1};
      
      % deal with partially observed observations
      if any(isnan(X(:)))
        
        C = obj.params.C;
        R = obj.params.R;
        
        G5 = zeros(M,M);
        for t=1:T
          
          XX = X(:,t) * X(:,t)';
          nidx = isnan(XX(:));
          if any(nidx)
            E = X(:,t) * Y(:,t)' * C' + R;
            XX(nidx) = E(nidx);
          end
          nidx = isnan(XX(:));
          if any(nidx)
            E = C * Y(:,t) * X(:,t)' + R;
            XX(nidx) = E(nidx);
          end
          nidx = isnan(XX(:));
          if any(nidx)
            E = C * (V{t} + Y(:,t)*Y(:,t)') * C' + R;
            XX(nidx) = E(nidx);
          end
          
          G5 =  G5 + XX;
          
        end
        
        G6 = zeros(M,K);
        for t=1:T
          
          XY = X(:,t) * Y(:,t)';
          EXY = C * (V{t} + Y(:,t)*Y(:,t)');
          nidx = isnan(XY(:));
          XY(nidx) = EXY(nidx);
          G6 = G6 + XY;
          
        end
        
      else
        G5 = X * X';
        G6 = X * Y';
      end
      
    end
    
    function LL = compute_loglik(obj,X,Y,V)
  
      T = size(X,2);

      C = obj.params.C;
      R = obj.params.R;
      
      LL = 0;
      for t=1:T
        
        obs = ~isnan(X(:,t));
        
        e = X(obs,t) - C(obs,:)*Y(:,t); % innovation (measurement error)
        S = C(obs,:)*V{t}*C(obs,:)' + R(obs,obs); % covariance of the innovation
        LL = LL - log(det(S)) ./ 2 - ((e' / S) * e) ./ 2;      
        LL = LL - (numel(obs)/2) * log(2*pi);
    
      end
      
      
    end
    
  end
  
end
