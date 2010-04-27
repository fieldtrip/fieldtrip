classdef rls < regressor
% Recursive Least Squares estimator
%
% Calculates regularised least squares regression in an online fashion,
% that is, it minimises  E=(Y-X*beta) + L2*|beta|^2 given one row of X 
% and Y at a time. Batch solution is  beta = (X'*X + L2)\Y of course, and
% RLS works by rank-1 updates to the inverse of (X'*X + L2).
%
% You can also set a forgetting factor 'lambda' (which is actually a 
% keep-in-mind-factor), which down-weights old samples. This also downweights
% the ridge parameter L2 over time (but that's fine if you have lots of data).
%
% Copyright (c) 2010, S.Klanke

  properties
	% invH and beta are part of 'params'
	% invH   = [];	% Inverse of regularised hessian   (X'*X
	% beta   = [];
	L2     = 1e-6;
    lambda = 1;
  end
    
    methods
        
        function obj = rls(varargin)
            obj = obj@regressor(varargin{:});
			if obj.lambda<=0 | obj.lambda>1
				error 'Forgetting factor lambda out of range';
			end
			if obj.L2 <=0
				error 'Regularisation (ridge) parameter must be positive';
			end
        end        
        
        function p = estimate(obj,X,Y)
			[N,dx]=size(X);
			[NY,dy]=size(Y);
			if N~=NY
				error 'X and Y must have the same number of rows';
			end

			if isempty(obj.params)
				invH = eye(dx) / obj.L2;
				beta = zeros(dx,dy);
			else
				invH = obj.params.invH;
				beta = obj.params.beta;
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
			p.invH = invH;
			p.beta = beta;
		end
        
        function Y = map(obj,X)
			Y = X*obj.params.beta;
        end
 
        function [m,desc] = getmodel(obj)
			% SK: not sure whether this matches the desired semantics
			m = {obj.params.invH obj.params.beta}
			desc = {'Inverse Hessian' 'beta'};
        end
    end
end 
