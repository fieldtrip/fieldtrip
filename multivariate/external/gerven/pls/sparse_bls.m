function [A,B,C,Ypredict] = sparse_bls(X,Y,nhidden,lambda1,lambda2,VERBOSE)

% SPARSE_BLS  Sparse orthogonalized partial least squares

% X: ninput x nsamples input data matrix
% Y: noutput x nsamples output data matrix
% nhidden: number of components [1]
% lambda1: l1 regularization parameter [1]
% lambda2: l2 regularization parameter [0.01]
%
% A: noutput x nhidden weight matrix
% B: ninput x nhidden weight matrix
% C: 1 x nhidden bias vector
% Ypredict: noutput x nsamples predictions


% Parse inputs

if nargin < 3,
    nhidden = 1;
end

if nargin < 4,
	lambda1 = 1;
end

if nargin < 5,
	lambda2 = 0.01;
end

ninput = size(X,1);
noutput = size(Y,1);
nsamples = size(X,2);
options = struct('offset',1,'maxiter',1e4,'tol',1e-2);

if nhidden > 1,      
    % Run nhidden times sequentially
    
    R = Y;
    A = zeros(noutput,nhidden);
    B = zeros(ninput,nhidden);
    C = zeros(1,nhidden);
    for i=1:nhidden,
        [A(:,i),B(:,i),C(i)] = sparse_bls(X,R,1,lambda1,lambda2,VERBOSE);
        if i < nhidden,    % deflate
            Z = B(:,i)'*X + C(i); % hidden activations
            R = R - A(:,i)*Z; % Y activations after deflation
        end
        if VERBOSE,
            fprintf('done %d out of %d; sparsity %g\n',i,nhidden,length(find(B(:,i)))/ninput);
        end
    end
else

    % Run for a single hidden unit
     
    B = zeros(ninput,1);
    C = 0;
    iter = 0;
    maxiter = 1000;
    tol = nhidden*noutput*(1e-10);
   
    % Initialize A to first principal component of Y
   
    optseig.disp = 0;
    if nsamples < noutput,
        [d1,d2] = eigs(Y'*Y,[],1,'LM',optseig);
        A = Y*d1;
        A = A/sqrt(A'*A);
    else
        [A,d2] = eigs(Y*Y',[],1,'LM',optseig);
    end
   
%   A = randn(size(A));   % random
%   A = A/sqrt(A'*A);

    Aold = A;
    while iter < maxiter,
	   
        if VERBOSE > 1,
            fprintf('   now starting iteration %d\n',iter+1);
        end
       
        Z = A'*Y;    % reconstruct Z given A from output Y
       
        % Use elastic net code to fit reconstructed Z given input X
        
        [B,C] = elastic(X,Z,lambda1,lambda2,options,B,C);

%         opts = glmnetSet;
%         opt.alpha = 1;
%         opts.lambda = lambda1;
%         res = glmnet(X',Z','gaussian',opts);
%         B = res.beta;
%         C = res.a0;
%         m = ft_mv_glmnet('validator',ft_mv_crossvalidator('nfolds',0.66,'metric','coefdet'),'alpha',1,'family','gaussian');
%         m = m.train(X',Z');
%         B = m.weights(1:(end-1));
%         C = m.weights(end);
        
        Z = B'*X + C;   % reconstruct Z given B and C
        
        % Find optimal A under constraint A'*A = 1
        
        Syz = Y*Z'/nsamples;
        denom = sqrt(Syz'*Syz);
        if denom,
            A = Syz/denom;
        else
            A = Aold;
        end
        if any(isnan(A(:))),    % check...
            error('nans!!!\n');
        end
           
        if sumsqr(A - Aold) < tol,
            if VERBOSE > 1,
                fprintf('   done!!\n',iter+1);
            end
		    iter = maxiter;
        else
	        iter = iter + 1;
            Aold = A;
        end
    end
end


% Parse outputs

if nargout > 3,
    Z = B'*X + C'*ones(1,nsamples);
    Ypredict = A*Z;
end