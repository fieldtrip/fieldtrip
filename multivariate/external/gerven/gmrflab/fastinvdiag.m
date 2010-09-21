function S = fastinvdiag(L)

% L: cholesky factorization of precision matrix
%
% S: diagonal of covariance matrix

N = length(L);

S = zeros(N,1);
tic
for k=1:N,
   e = spalloc(N,1,1);
   e(k) = 1;
   x = L\e;
   S(k) = x'*x;
   if ~rem(k,1000),
      fprintf('%d out of %d (%g %%); %g\n',k,N,k/N*100,toc);
      tic
   end
end