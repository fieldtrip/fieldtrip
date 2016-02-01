function [A,B,C,G] = ospls(X,Y,nhidden,opts)

% OSPLS  Sparse orthogonalized partial least squares where the elastic
% net is implemented using either a native implementation or the glmnet package
%
% The native implementation supports arbitrary alpha (ridge penalty) and
% lambda (lasso penalty) values and supports a coupling matrix for alpha to
% implement smoothing. It requires a fixed value for lambda.
%
% The glmnet implementation is much faster but uses alpha to mix between
% ridge and lasso while lambda determines that amount of regularization. It
% also does not support a coupling matrix. It learns the optimal lambda
% using inner cross-validation.
%
% X: ninput x nsamples input data matrix
% Y: noutput x nsamples output data matrix
% nhidden: number of components [1]
% opts: glmnet parameters
%
% A: noutput x nhidden weight matrix
% B: ninput x nhidden weight matrix
% C: 1 x nhidden bias vector
% G: output parameters for debugging purposes


% Parse inputs

if nargin < 3,
    nhidden = 1;
end

ninput = size(X,1);
noutput = size(Y,1);
nsamples = size(X,2);

if nhidden > 1,     
  
    G = cell(1,nhidden);
  
    % Run nhidden times sequentially
    
    R = Y;
    A = zeros(noutput,nhidden);
    B = zeros(ninput,nhidden);
    C = zeros(1,nhidden);
    for i=1:nhidden,
        [A(:,i),B(:,i),C(i),G{i}] = ospls(X,R,1,opts);
        if i < nhidden,    % deflate
            Z = B(:,i)'*X + C(i); % hidden activations
            R = R - A(:,i)*Z; % Y activations after deflation
        end
        if opts.verbose,
            fprintf('done %d out of %d; proportion selected %g\n',i,nhidden,nnz(B(:,i))/ninput);
        end
    end
else

    % Run for a single hidden unit
     
    B = zeros(ninput,1);
    C = 0;
    iter = 0;
    maxiter = 100;
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
      
    glmnetopts = glmnetSet;
    
    g = ft_mv_glmnet('method',opts.method,'family','gaussian','validator',ft_mv_crossvalidator('nfolds',5,'metric','correlation'));
    
    FN = fieldnames(glmnetopts);
    for c=1:length(FN)
      g.(FN{c}) = glmnetopts.(FN{c});
    end
    FN = fieldnames(opts);
    for c=1:length(FN)
      if isfield(glmnetopts,FN{c})
        g.(FN{c}) = opts.(FN{c});
      end
    end
    
    Aold = A;
    while iter < maxiter,
	   
        if opts.verbose > 1,
            fprintf('   now starting iteration %d\n',iter+1);
        end
       
        Z = A'*Y;    % reconstruct Z given A from output Y
       
        % Use glmnet elastic net code to fit reconstructed Z given input X
        try
        
        
          f = g.train(X',Z');
          
          B = f.weights(1:(end-1));
          C = f.weights(end);
          
        catch
          
          % we end up in the catch block when the number of included
          % variables pmax has been exceeded; this should not happen
          warning(lasterr);
          
          B = zeros(ninput,1);
          C = 0;
          
          % this implies that deflation has no effect anymore
          
        end
        
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
            if opts.verbose > 1,
                fprintf('   done!!\n',iter+1);
            end
		    iter = maxiter;
        else
	        iter = iter + 1;
            Aold = A;
        end
    end

    G = f;
    
end
