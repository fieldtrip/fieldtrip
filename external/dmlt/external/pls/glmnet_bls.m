function [A,B,C,G] = glmnet_bls(X,Y,nhidden,opts,VERBOSE)

% GlMNET_BLS  Sparse orthogonalized partial least squares where the elastic
% net is implemented using the glmnet package
%
% X: ninput x nsamples input data matrix
% Y: noutput x nsamples output data matrix
% nhidden: number of components [1]
% opts: glmnet parameters
%
% A: noutput x nhidden weight matrix
% B: ninput x nhidden weight matrix
% C: 1 x nhidden bias vector
% G: glmnet objects for debugging purposes


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
        [A(:,i),B(:,i),C(i),G{i}] = glmnet_bls(X,R,1,opts,VERBOSE);
        if i < nhidden,    % deflate
            Z = B(:,i)'*X + C(i); % hidden activations
            R = R - A(:,i)*Z; % Y activations after deflation
        end
        if VERBOSE,
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
   
    Aold = A;
    while iter < maxiter,
	   
        if VERBOSE > 1,
            fprintf('   now starting iteration %d\n',iter+1);
        end
       
        Z = A'*Y;    % reconstruct Z given A from output Y
       
        % Use glmnet elastic net code to fit reconstructed Z given input X
        try
        
          f = ft_mv_glmnet('family','gaussian','validator',ft_mv_crossvalidator('nfolds',5,'metric','correlation'));
          FN = fieldnames(opts);
          for c=1:length(FN)
            f.(FN{c}) = opts.(FN{c});
          end
          
          f = f.train(X',Z');
          
          B = f.weights(1:(end-1));
          C = f.weights(end);
        
        catch
          
          % we end up in the catch block when the number of included
          % variables pmax has been exceeded; this should not happen
          warning('maximal number of variables exceeded');
          
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
            if VERBOSE > 1,
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
