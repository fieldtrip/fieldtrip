classdef slda < dml.method
% SLDA shrinkage linear discriminant analysis.
%
%   DESCRIPTION
%   Uses shrinkage to estimate covariance matrices:
%   'diaguneq' : Target D of Schafer and Strimmer; Software of Strimmer,
%                covshrinkKPM(X(Y==k,:),0) would give the same result
%                (default)
%   'diagcommon' : Target B of Schafer and Strimmer as in Blankertz et al 
%
%   REFERENCE
%   J. Schafer and K. Strimmer (2005) A Shrinkage Approach to Large-Scale
%   Covariance Matrix Estimation and Implications for Functional Genomics
%  
%   Blankertz B, Lemm S, Treder M, Haufe S, Müller K-R. Single-Trial
%   Analysis and Classification of ERP Components - a Tutorial. Neuroimage.
%   2010
%
%   Opgen-Rhein R, Strimmer K. Accurate ranking of differentially
%   expressed genes by a distribution-free shrinkage approach. Stat Appl Genet Mol Biol. 2007;6:Article9. 
%
%   EXAMPLE
%
%   NOTE
%   The approach by Opgen-Rhein may also be of use (see covshrinkKPM)
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)

  properties
    
    pi % priors
    
    mu % means
    
    Sigma % covariances
    
    Sinv % inverse of joint covariance matrix
    
    lambda % regularization parameter (automatically determined)
    
    shrinkage = 'diaguneq'; % 'diaguneq' or 'diagcommon'
    
  end
  
  methods
    
    function obj = slda(varargin)

      obj = obj@dml.method(varargin{:});

    end
    
    function obj = train(obj,X,Y)
      
      nclasses = numel(unique(Y));
      d = size(X,2);
      
      % estimate class priors
      obj.pi = zeros(nclasses,1);
      for j=1:nclasses
        obj.pi(j) = sum(Y==j)/size(Y,1);
      end

      % estimate class-conditional means
      obj.mu = cell(nclasses,1);
      for k=1:nclasses
        obj.mu{k} = nanmean(X(Y == k,:))';
      end
         
      % estimate class-conditional covariance matrices
      obj.Sigma = cell(nclasses,1); 
      if isempty(obj.lambda)
        obj.lambda = zeros(nclasses,1);
      end
      for k=1:nclasses
        
        obj.Sigma{k} = cov(X(Y == k,:));
        
        if obj.lambda(k)==0
          
          mX = bsxfun(@minus,X(Y == k,:),obj.mu{k}');
          N = size(mX,1);
          W = zeros(N,d,d);
          for n=1:N
            W(n,:,:) = mX(n,:)' * mX(n,:);
          end
          WM = mean(W,1);
          S = squeeze((N/(N-1)) .* WM);
          
          switch obj.shrinkage
            
            case 'diagcommon'
              % Target 'B' of Schafer and Strimmer; choice by Blankertz et al.
              
              VS = squeeze((N/((N-1).^3)) .* sum(bsxfun(@minus,W,WM).^2,1));
              
              v = mean(diag(S));
              
              t = triu(S,1);
              obj.lambda(k) = sum(VS(:)) / (2*sum(t(:).^2) + sum((diag(S)-v).^2));
              
            case 'diaguneq'
              % Target 'D' of Schafer and Strimmer
              
              mX = zscore(X(Y==k,:));
              N = size(mX,1);
              W = zeros(N,d,d);
              for n=1:N
                W(n,:,:) = mX(n,:)' * mX(n,:);
              end
              WM = mean(W,1);
              VS = squeeze((N/((N-1).^3)) .* sum(bsxfun(@minus,W,WM).^2,1));
              
              t = triu(VS,1); u = triu(corr(mX),1);
              obj.lambda(k) = sum(t(:)) / sum(u(:).^2);
              
          end
          
          obj.lambda(k) = max(0,min(1,obj.lambda(k)));
          
        end
        
        % the regularizing matrix
        switch obj.shrinkage
          
          case 'diagcommon'
            % Target 'B' of Schafer and Strimmer; choice by Blankertz et al.
            
            T = v*eye(d);
            
          case 'diaguneq'
            % Target 'D' of Schafer and Strimmer
            
            T = diag(diag(obj.Sigma{k}));
            
        end
        
        obj.Sigma{k} = (1-obj.lambda(k))*obj.Sigma{k} + obj.lambda(k)*T;

      end
      
      % compute inverse of joint covariance matrix (i.e. LDA instead of QDA)
      S = 0; for c=1:length(obj.Sigma), S = S + obj.Sigma{c}; end
      obj.Sinv = inv(S/length(obj.Sigma));
      
    end
        
    
    function Y = test(obj,X)
      
      % compute discriminant functions d
      d = zeros(size(X,1),length(obj.Sigma));
      for k=1:length(obj.Sigma)
        d(:,k) = X * obj.Sinv * obj.mu{k} - obj.mu{k}' * obj.Sinv * obj.mu{k} + log(obj.pi(k));
      end
      
      % get most likely class
      Y = (d == repmat(max(d,[],2),[1 size(d,2)]));
      Y = bsxfun(@rdivide,Y,sum(Y,2));
      
    end

  end
  
end
