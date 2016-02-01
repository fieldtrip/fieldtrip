function S = hammersley(D,N)
% HAMMERSLEY - Hammersley quasi-random sequence
% 
%   S = HAMMERSLEY(D,N) Generates N numbers from a D-dimensional
%   Hammersley quasirandom sequence using successive primes
%   as bases except for the last dimension which has the form
%   [1:N]'/N-1/2/N (where the last term modifies usual Hammersley
%   sequence to produce sequence in open interval (0,1)). The 
%   matrix S is of size DxN.
  
% Copyright (c) 2008 Aki Vehtari

% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

S=zeros(N,D);
% Last column is easy
S(:,D)=([1:N]'/N)-1/(2*N);
% Matlab's PRIMES returns primes smaller than given value,
% but we need D-1 primes
pn=2*D;
p=primes(pn);
while(numel(p)<D-1)
  pn=2*pn;
  p=primes(pn);
end
P=p(1:D-1);
for k=1:D-1 % dimensions
  pk=P(k);
  for j=1:N     % points in sequence
    bj=j;
    n = max(1,round(log2(bj+1)/log2(pk)));
    while pk.^n <= bj
      n = n + 1;
    end
    b=zeros(1,n);
    b(n)=rem(bj,pk);
    while bj && n>1
      n=n-1;
      bj=floor(bj/pk);
      b(n)=rem(bj,pk);
    end
    S(j,k)=sum(fliplr(b)./pk.^[1:numel(b)]);
  end
end

S = S';