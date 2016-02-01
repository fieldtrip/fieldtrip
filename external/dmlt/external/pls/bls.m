function [A,B,Ypredict] = bls(X,Y,nhidden,algorithm)
% BLS  Bottleneck least squares aka sparse orthogonalized least squares
%
% X: ninput x nsamples input data matrix
% Y: noutput x nsamples output data matrix
% nhidden: dimension of hidden
%
% A: noutput x nhidden weight matrix
% B: ninput x nhidden weight matrix


if nargin < 3,
  nhidden = min(size(Y,1),2);
end

if nargin < 4,
  algorithm = 'direct';
end

[ninput,nsamples] = size(X);

Sxx = symmetric(X*X')/nsamples;

switch algorithm
  
  case 'direct'
    Syx = Y*X'/nsamples;
    Sxy = Syx';
    Sxxinv = symmetric(pinv(Sxx));
    S = symmetric(Syx * Sxxinv * Sxy);
    opts.disp = 0;
    opts.issym = 'true';
    [A,D] = eigs(S,nhidden,'LM',opts);   % first nhidden eigenvectors of S
    B = Sxxinv*Sxy*A;
    
    
  case 'indirect'          % same result and subspace as 'direct', but then different linear combination
    
    Syx = Y*X'/nsamples;
    Sxy = Syx';
    Sxxinv = symmetric(pinv(Sxx));
    
    % initialize randomly
    
    B = randn(ninput,nhidden);
    
    iter = 0;
    maxiter = 1000;
    tol = nhidden*ninput*(1e-10);
    Bold = B;
    while iter < maxiter,
      Z = B'*X;
      %		Szz = symmetric(Z*Z')/nsamples;
      Szz = symmetric(B'*Sxx*B);
      %		Syz = symmetric(Y*Z')/nsamples;
      Syz = Syx*B;
      A = Syz*pinv(Szz);
      
      B = Sxxinv*Sxy*A*pinv(A'*A);
      if sumsqr(B - Bold) < tol,
        fprintf('done after %d iterations!!\n',iter+1);
        iter = maxiter;
      else
        iter = iter + 1;
        Bold = B;
      end
    end
    
  case 'sequential'      % same result as 'direct', up to minus sign per hidden unit
    
    R = Y;
    noutput = size(Y,1);
    A = zeros(noutput,nhidden);
    B = zeros(ninput,nhidden);
    for i=1:nhidden,
      [A(:,i),B(:,i)] = bls(X,R,1,'direct');
      if i < nhidden,
        R = R - A(:,i)*B(:,i)'*X;
      end
    end
    
end

% transform to make latent variable orthonormal

Szz = symmetric(B'*Sxx*B);
Szzinv = symmetric(pinv(Szz));
sqrtSzzinv = chol(Szzinv(end:-1:1,end:-1:1));    % such that the rescaled B(i,:) is a linear combi of the unscaled B(1:i,:); this only makes sense for the direct method
sqrtSzzinv = sqrtSzzinv(end:-1:1,end:-1:1);
B = B*sqrtSzzinv;
A = A*pinv(sqrtSzzinv);

if nargout > 2,
  Ypredict = A*B'*X;
end

function A = symmetric(A)

A = (A+A')/2;
	