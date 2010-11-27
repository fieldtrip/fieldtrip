classdef ft_mv_rls < ft_mv_predictor
% Recursive Least Squares estimator
%
% Calculates regularised least squares regression in an online fashion,
% that is, it minimises  E=(Y-X*beta) + L2*|beta|^2 given one row of X 
% and Y at a time. Batch solution is  beta = inv(X'*X + L2)*X'*Y of course, and
% RLS works by rank-1 updates to the inverse of (X'*X + L2) and similarly for
% 'beta' itself.
%
% You can set a forgetting factor 'lambda' (actually a keep-in-mind-factor), 
% which down-weights old samples. This also downweights the ridge parameter 
% L2 over time (but that's fine if you have lots of data).
%
% Copyright (c) 2010, S.Klanke

  properties
	  
    L2     = 1e-6;
    lambda = 1;
  
    invH   = [];	% Inverse of regularised hessian   (X'*X + L2*eye)
    beta   = [];
    
  end
    
  methods
    
    function obj = ft_mv_rls(varargin)
      obj = obj@ft_mv_predictor(varargin{:});
      if obj.lambda<=0 || obj.lambda>1
        error 'Forgetting factor lambda out of range';
      end
      if obj.L2 <=0
        error 'Regularisation (ridge) parameter must be positive';
      end
    end
    
    function obj = train(obj,X,Y)
      
      [N,dx]=size(X);
      [NY,dy]=size(Y);
      if N~=NY
        error 'X and Y must have the same number of rows';
      end
      
      if isempty(obj.invH)
%        invH = eye(dx) / obj.L2;
%        beta = zeros(dx,dy);
        obj.invH = inv(X'*X + obj.L2*eye(size(X,2),size(X,2)));
        obj.beta = inv(X'*X + obj.L2*eye(size(X,2),size(X,2)))*X'*Y;
        return;
      else
        invH = obj.invH;
        beta = obj.beta;
      end
      
      lam = obj.lambda;
      for n=1:N
        x = X(n,:); % keep as row vector
        y = Y(n,:); % keep as row vector
        invHx = invH*x';
        L = invHx/(lam+x*invHx);
        invH = invH - L*invHx';
        beta = beta + L*(y - x*beta);
      end
      obj.invH = invH;
      obj.beta = beta;
      
    end
    
    function Y = test(obj,X)
      Y = X*obj.beta;
    end
    
    function [m,desc] = model(obj)
      m = {obj.invH obj.beta};
      desc = {'Inverse Hessian' 'beta'};
    end
    
  end
end
