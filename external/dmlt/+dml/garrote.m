classdef garrote < dml.method
% GARROTE variational garrote.
%
%   DESCRIPTION
%   Variational garrote implementation of sparse regression
%
%   REFERENCE
%   http://arxiv.org/abs/1109.0486
%
%   EXAMPLE
%   n=10; % input dimension 
%   p=100; % number of samples
%   pt = 2*n; 
%   ns = round(.05*n);  % sparsity level
%   betax=1;            % inverse noise in the input
%   betah=1;            % inverse noise response variance
%   snonzero=randperm(n);
%   snonzero(ns+1:end) = [];
%   w=sparse(1,n);
%   w(snonzero)=1;
%   sigma=sqrt(1/betah);    % noise response
%   sigmax=sqrt(1/betax);    % noise input  
%   x=sigmax*randn(n,p);
%   x=x-mean(x,2)*ones(1,p);
%   dx=sqrt(1/p*sum(x.^2,2));
%   x=x./(dx*ones(1,p));
%   y=w*x+sigma*randn(1,p);
%   y=y-mean(y);
%   xt=sigmax*randn(n,pt);
%   xt=xt-mean(xt,2)*ones(1,pt);
%   yt=w*xt+sigma*randn(1,pt);
%   yt=yt-mean(yt);
%   m = dml.garrote('valset',20);
%   m = m.train(x',y');
%   yp = m.test(xt');
%   plot([yt; yp]');
%
%   DEVELOPERS
%   Bert Kappen (b.kappen@science.ru.nl)
%   Vicenc Gomez (v.gomez@science.ru.nl)


  properties
    
     method='regression'  % method for optimization 'dual' or 'regression' for fixed gamma ('dual')
     maxiter=1e4    % maximum number of iterations for optimization for fixed gamma (1e4)
     max_sum_m      % increases gamma values until sum(m)=max_sum_m  (n/2)
     beta_max=1e3   % increases gamma values until beta=beta_max (1e3)
     n_gamma=50     % number of gamma values to scan (50)
     dmmin=1e-12    % convergence threshold for mean field error (1e-12)
     valset         % part of the training set used for validation (0.1*p)
    
     res            % learnt parameters
     
  end
  
  methods
    
    function obj = garrote(varargin)

      obj = obj@dml.method(varargin{:});

    end
    
    function obj = train(obj,X,Y)
        
      [n,p] = size(X);
      if isempty(obj.valset)
        obj.valset=0.1*n;
      end
      
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
      
      if isempty(obj.max_sum_m), obj.max_sum_m = n/2; end
      
      warning off
      args=struct(obj);    
      warning on
      obj.res = vg_method_train(X',Y',args);

    end
    
    function Y = test(obj,X)
     % return predicted Y based on input X
     
     Y = vg_method_test(X',obj.res);
      
    end

    function m = model(obj)

      m.m_mf = obj.m_mf;
      
    end

  end
  
end
