function [itr, itst] = cvit(n, k, rsubstream)
%CVIT  Create itr and itst indeces for k-fold-cv
%
%  Description
%    [ITR,ITST]=CVIT(N,K) returns 1xK cell arrays ITR and ITST
%    holding cross-validation indeces for train and test sets
%    respectively. K-fold division is balanced with all sets having
%    floor(N/K) or ceil(N/K) elements.
%
%    [ITR,ITST]=CVIT(N,K,RS) with integer RS>0 also makes random
%    permutation, using substream RS. This way different
%    permutations can be produced with different RS values, but
%    same permutation is obtained when called again with same RS. 
%    Function restores the previous random stream before exiting.
%

% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

if nargin<3
  rsubstream=0;
end
if nargin < 2
  k=10;
end
a=k-rem(n,k);
b=floor(n/k);
for cvi=1:a
  itst{cvi}=[[1:b]+(cvi-1)*b]; 
  itr{cvi}=setdiff(1:n,itst{cvi}); 
end
for cvi=(a+1):k
  itst{cvi}=[(a*b)+[1:(b+1)]+(cvi-a-1)*(b+1)]; 
  itr{cvi}=setdiff(1:n,itst{cvi}); 
end  

if rsubstream>0
  prevstream=setrandstream(0, 'mrg32k3a');
  stream.Substream = rsubstream;

  rii=randperm(n);
  for cvi=1:k
    itst{cvi}=rii(itst{cvi});
    itr{cvi}=rii(itr{cvi});
  end

  setrandstream(prevstream);
end
