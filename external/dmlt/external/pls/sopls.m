function [A,B,C,G] = sopls(X,Y,nhidden,reg,verbose)
% SOPLS  Sparse orthogonalized partial least squares where the elastic
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
% reg: employed regression method
%
% A: noutput x nhidden weight matrix
% B: ninput x nhidden weight matrix
% C: 1 x nhidden bias vector
% G: output parameters for debugging purposes


% Parse inputs

if nargin < 3, nhidden = 1; end
if nargin < 4, error('please specify regularizer'); end
if nargin < 5, verbose = false; end

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
      
      if verbose
        fprintf('learning hidden variable %d of %d\n',i,nhidden);
      end
      
      [A(:,i),B(:,i),C(i),G{i}] = sopls(X,R,1,reg);
      if i < nhidden,    % deflate
        Z = B(:,i)'*X + C(i); % hidden activations
        R = R - A(:,i)*Z; % Y activations after deflation
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
	          
        Z = A'*Y;    % reconstruct Z given A from output Y
 
        try

          f = reg.train(X',Z');
          
          % all regularizers should use this convention
          B = f.model.weights;
          C = f.model.bias;
                    
        catch
          
          % this should not happen
          warning(lasterr);
          B = zeros(ninput,1);
          C = 0;
          
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
          iter = maxiter;
        else
          iter = iter + 1;
          Aold = A;
        end
    end

    G = f;
    
end
