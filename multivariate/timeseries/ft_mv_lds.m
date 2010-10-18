classdef ft_mv_lds < ft_mv_timeseries
%FT_MV_LDS linear dynamical system 
%
% state can be partially observed/unobserved during training.
% partial observability of multiple observations is supported
% observations are normally distributed conditional on the state.
% multiple observations sequences are supported
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
%  epsilon = 0 % added to the diagonal of the covariance matrices for numerical
%                stability (e.g. 1e-7);
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
%
% NOTES: 
% - X and Y are swapped wrt the Kalman filter conventions
% - obj.A = zeros(K,K) will remove the dynamics from the model
%
% EXAMPLE:
%
% rand('seed',3); randn('seed',3);
% 
% nsamples = 1000; ncov = 10; ncycles = 10;
% Y = sin(ncycles * 2 * pi * (1:nsamples) ./ nsamples)';
% 
% k = ft_mv_lds('inference','smooth','verbose',true);
% X = repmat(Y,[1 ncov]) + randn(size(Y,1),ncov); 
% k = k.train(zscore(X),zscore(Y)); % everything assumed observed
% X = repmat(Y,[1 ncov]) + randn(size(Y,1),ncov); 
% Z = k.test(zscore(X));
% plot(Z,'ko');
% disp(mean(abs(Z - Y)));
% 
% k = ft_mv_lds('verbose',true);
% k = k.train(zscore(X),nan(size(Y))); % hidden state estimation
% X = repmat(Y,[1 ncov]) + randn(size(Y,1),ncov); 
% Z = k.test(zscore(X));
% plot(Z,'ko');
% disp(mean(abs(Z - Y)));
% 
%
% Copyright (c) 2010, Marcel van Gerven, Ali Bahramisharif
  
  properties
    
    maxiter = 100; % maximum number of EM iterations
    thresh = 1e-4; % EM convergence threshold    
    diagQ = 0; % diagonal state noise
    diagR = 0; % diagonal measurement noise
    
    inference = 'smooth'; % filter / smooth
    
    epsilon = 0; % numerical stability
    
    A
    C
    mu0
    V0
    R
    Q
    
    loglik
    
  end

  methods
    
    function obj = ft_mv_lds(varargin)
      
      obj = obj@ft_mv_timeseries(varargin{:});
    end
    
    function obj = train(obj,X,Y)
               
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
      eidx = eidx & cellfun(@(x)(isempty(x)),Y);
      X = X(~eidx);
      Y = Y(~eidx);
        
      K = size(Y{1},1);
      M = size(X{1},1);
      N = length(X);
    
      if isempty(obj.A), obj.A = 1e-3*randn(K,K); end
      if isempty(obj.C), obj.C = 1e-3*randn(M,K); end
      if isempty(obj.mu0), obj.mu0 = 1e-3*randn(K,1); end
      if isempty(obj.V0), obj.V0 = 1e-3*eye(K); end
      if isempty(obj.R), obj.R = 1e-3*eye(M); end
      if isempty(obj.Q), obj.Q = 1e-3*eye(K); end
        
%       % DEBUG
%       [obj.A, obj.C, obj.Q, obj.R, obj.mu0, obj.V0, obj.loglik] = learn_kalman(X,obj.A,obj.C,obj.Q,obj.R,obj.mu0,obj.V0);
%       return;
      
      oldLL = inf;
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
          fprintf('EM step: %d; log likelihood: %g\n',iter,LL);
        end
        
        if LL < oldLL, fprintf('non-decreasing log likelihood!\n'); end
  
        % M step
        T = sum(cellfun(@(x)(size(x,2)),X));
        
        obj.C = G6 / G1;
        
        obj.R = (G5 - obj.C * G6') ./ T;
        obj.R = (obj.R + obj.R') ./ 2;
        
        if obj.diagR
          obj.R = (1-obj.diagR) * obj.R + obj.diagR * diag(diag(obj.R));
        end
        
        if ~all(obj.A==0) 
          obj.A = G4 / G3;
        end
        
        obj.Q = (G2 - obj.A * G4') ./ (T-N);
        obj.Q = (obj.Q + obj.Q') ./ 2;
        
        if obj.diagQ
          obj.Q = (1-obj.diagQ) * obj.Q + obj.diagQ * diag(diag(obj.Q));
        end
        
        obj.mu0 = mu0/N;
        obj.V0 = V0/N - obj.mu0*obj.mu0';      
        
        % symmetricize and make positive semidefinite
        obj.V0 = (obj.V0 + obj.V0') ./ 2;
        obj.V0 = obj.V0 + obj.epsilon*eye(size(obj.V0));
         
        iter = iter + 1;
        
        if ~any(isnan(Y{1}(:)))
          break;
        end
        
      end
      
      if obj.verbose
        fprintf('EM step: %d; loglikelihood: %g\n',iter,LL);
      end
      
      obj.loglik = loglik(2:end);
          
    end
    
    function [mu,V,loglik] = test(obj,X)   
     % LDS inference

     if ~iscell(X), X = {X}; end
     
     nsets = length(X);
      
     mu = cell(1,nsets);
     for c=1:nsets
     
       % representation as variables * timepoints
       XX = X{c}';
       
       [m,V,loglik] = obj.filter(XX);
       
       if strcmp(obj.inference,'smooth')
         [m,V] = obj.smooth(m,V);
       end
       
       mu{c} = m';
     
     end
     
     if length(nsets)==1
       mu = mu{1};
     end
     
    end
    
    function [mu,V,J,VP] = smooth(obj,mu1,V1,Y)
      % Kalman smoother; uses filtered means and variances
      
      T = numel(V1);

      nandim=[];
      dim=[];
      if exist('Y','var')
          for i=1:size(Y,1)
              if isnan(Y(i,1))
                  nandim=[nandim i];
              else
                  dim=[dim,i];
              end
          end

          A = obj.A(nandim,nandim);
          C = obj.C(:,nandim);
          Q = obj.Q(nandim,nandim);
      else

          A   = obj.A;
          Q   = obj.Q;
          C   = obj.C;
      end
      % P(t,t-1) at horizon
      P = A * V1{T-1} * A' + Q;
      PC = P * C';
      K = PC / (C * PC + obj.R);
      VP=cell(1,T);
      KC=K*C;
      VP{T}=(eye(size(KC))-KC)*A*V1{T-1};
      
      % RTS equations
      
      mu = mu1;
      V = V1;
      J = cell(1,T-1); % needed to compute transition matrix A in EM
      for n=T:-1:2

        P = A * V1{n-1} * A' + Q;
        
        if det(P)<0
          disp('bad predictive covariance')
          lambda=max(svd(P));
          P=P+(lambda+obj.epsilon)*eye(size(P));
        end
        
        if ~all(V{n}(:)==0)
          
          J{n-1} = (V1{n-1} * A') / P;
          
          mu(:,n-1) = mu1(:,n-1) + J{n-1} * (mu(:,n) - A * mu1(:,n-1));
          V{n-1} = V1{n-1} + J{n-1} * (V{n} - P)*J{n-1}';

          if n<T
            VP{n}=V1{n}*J{n-1}'+J{n}*(VP{n+1}-A*V1{n})*J{n-1}';
          end
          
        else
          
          % fix when variance is zero (deterministic state)
          J{n} = zeros(size(V{n}));
          mu(:,n) = mu1(:,n);
          V{n} = zeros(size(V{n}));
          VP{n} = zeros(size(V{n}));
          
        end
        
        % symmetricize and make positive semidefinite
        V{n} = (V{n} + V{n}') ./ 2;
        V{n} = V{n} + obj.epsilon*eye(size(V{n}));
        
      end
      
    end
    
    function [mu,V,LL] = filter(obj,X,varargin)   
      % Kalman filter algorithm; returns mean mu and covariance matrix V of
      % the states, the log likelihood and the observations Xobs with missing
      % values imputed (the latter remains to be done)
      
      nandim=[];
      dim=[];
      if ~isempty(varargin)
          Y=varargin{1};
      end
      if exist('Y','var')
          for i=1:size(Y,1)
              if isnan(Y(i,1))
                  nandim=[nandim i];
              else
                  dim=[dim,i];
              end
          end
      
          A = obj.A(nandim,nandim);
          C = obj.C(:,nandim);
          mu0 = obj.mu0(nandim);
          V0 = obj.V0(nandim,nandim);
          R = obj.R;
          Q = obj.Q(nandim,nandim);
          X1=X;
          X=X-obj.C(:,dim)*Y(dim,:);
      else
          A = obj.A;
          C = obj.C;
          mu0 = obj.mu0;
          V0 = obj.V0;
          R = obj.R;
          Q = obj.Q;
          Y=zeros(size(A,1),size(X,2));
          X1=X;
      end          
      
      K = size(A,1);
      T = size(X,2);
      
      mu = nan(K,T);
      V = cell(1,T);
      I = eye(K);      

      Xpred=mu;
      Vpred=V; 

      % take measurements into account
      if isempty(X) || all(isnan(X(:,1)))
       
        mu(:,1) = mu0;
        V{1} = V0;
      
      else
        Xpred(:,1)=mu0;
        Vpred{1}=V0;
        obs = ~isnan(X(:,1));
        Cb = C(obs,:); % emission matrix for observed measurements
      
        % Kalman gain matrix
        K = (V0 * Cb') / (Cb * V0 * Cb' + R(obs,obs));

        mu(:,1) = mu0 + K * (X(obs,1) - Cb * mu0);
        V{1}  = (I - K * Cb) * V0;
      
      end
      
      
      % symmetricize and make positive semidefinite
      V{1} = (V{1} + V{1}') ./ 2;
      V{1} = V{1} + obj.epsilon*eye(size(V{1}));
     
      for t=2:T

        if numel(nandim)>0
          AM = A * mu(:,t-1)+obj.A(nandim,dim)*Y(dim,t-1);
        else
            AM = A * mu(:,t-1);
        end
        P = A * V{t-1} * A' + Q;
        
        if isempty(X) || all(isnan(X(:,t)))
          
          mu(:,t) = AM;
          V{t} = P;
        
        else
            Xpred(:,t)=AM;
            Vpred{t}=P;
          obs = ~isnan(X(:,t));
          Cb = C(obs,:); % emission matrix for observed measurements
          
          PC = P * Cb';
          
          K = PC / (Cb * PC + R(obs,obs));
          
          mu(:,t) = AM + K * (X(obs,t) - Cb * AM);
          V{t}  = (I - K * Cb) * P;
        
        end
        
        % symmetricize and make positive semidefinite
        V{t} = (V{t} + V{t}') ./ 2;
        V{t} = V{t} + obj.epsilon*eye(size(V{t}));
  
      end
      
      % compute log likelihood using error decomposition ignoring constant terms  
      if nargout >= 3
          V1 = cell(1,T);
          K=size(obj.A,1);
          Y1=zeros(K,T);
          Y1(dim,:)=Y(dim,:);
          if numel(nandim)<1
              Y1=mu;
              V1=V;
          else
              Y1(nandim,:)=Xpred;
              i=1;
              V1{i}=zeros(K);
              V1{i}(dim,dim)=obj.Q(dim,dim);
              V1{i}(dim,nandim)=obj.Q(dim,nandim);
              V1{i}(nandim,dim)=obj.Q(nandim,dim);
              V1{i}(nandim,nandim)=Vpred{i};
              for i=2:T
                  V1{i}=zeros(K);
                  V1{i}(dim,dim)=obj.Q(dim,dim);
                  V1{i}(dim,nandim)=obj.Q(dim,nandim);
                  V1{i}(nandim,dim)=obj.Q(nandim,dim);
                  V1{i}(nandim,nandim)=Vpred{i};
              end
          end
          LL = obj.compute_loglik(X1,Y1,V1);
      end
      
    end
        
  end
  
  methods(Access=protected)
    
    function [G1,G2,G3,G4,G5,G6,mu0,V0,LL] = Estep(obj,X,Y)
              
      K = size(Y,1);
      M = size(X,1);
      T = size(X,2);

      nandim=[];
      dim=[];
      if exist('Y','var')
          for i=1:K
              if isnan(Y(i,1))
                  nandim=[nandim i];
              else
                  dim=[dim,i];
              end
          end
      end
      
      if any(isnan(Y(:))) % hidden state

          [mu1,V1,LL] = obj.filter(X,Y);
          [Y1,V1,J1,VP1] = obj.smooth(mu1,V1,Y);
               
        
          
          Y(nandim,:)=Y1;
          VP=cell(1,T);
          V = cell(1,T);
          J = cell(1,T-1);
          for i=1:T
              if i>1
                  VP{i}=zeros(K);
                  VP{i}(nandim,nandim)=VP1{i};
              end
              V{i}=zeros(K);
              V{i}(dim,dim)=V{i}(dim,dim)+obj.epsilon*eye(length(dim));
              V{i}(nandim,nandim)=V1{i};
              if i<T
                  J{i}=zeros(K);
                  J{i}(nandim,nandim)=J1{i};
              end
          end


      else
      
        V = repmat({0},size(Y));
        J = repmat({0},[size(Y,1) size(Y,2)-1]);
        VP = repmat({0},size(Y));
        LL = nan;
        
      end
      


      
      G1 = zeros(K,K);
      G4 = zeros(K,K);
      for t=1:T
        
        G1  = G1 + V{t} + Y(:,t)*Y(:,t)';
        
        if t>1
        
          G4 = G4 + Y(:,t)*Y(:,t-1)' + VP{t};
          %G4 = G4 + J{t-1}*V{t} + Y(:,t)*Y(:,t-1)'; % bishop
        
        end
        
      end
      
      G2 = G1 - Y(:,1) * Y(:,1)' - V{1};
      G3 = G1 - Y(:,T) * Y(:,T)' - V{T};
      
      mu0 = Y(:,1);
      V0 = V{1};
      
      % deal with partially observed observations
      if any(isnan(X(:)))
        
        C = obj.C;
        R = obj.R;
        
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

      C = obj.C;
      R = obj.R;
      
      LL = 0;
      for t=1:T
        
        obs = ~isnan(X(:,t));
        
        e = X(obs,t) - C(obs,:)*Y(:,t); % innovation (measurement error)
        S = C(obs,:)*V{t}*C(obs,:)' + R(obs,obs); % covariance of the innovation
        LL = LL - log(det(S)) ./ 2 - ((e' / S) * e) ./ 2;
        if ~isreal(LL)
            LL=real(LL);
            disp('bad likelihood')
        end
        LL = LL - (numel(obs)/2) * log(2*pi); %length of the feature vector 
    
      end
      
      
    end
    
  end
  
end
