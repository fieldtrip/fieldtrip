function s = resampstr(p,m,n)
%RESAMPSTR Stratified resampling
%
%  Description
%    S = RESAMPSTR(P) returns a new set of indices according to 
%      the probabilities P. P is array of probabilities, which are
%      not necessarily normalized, though they must be
%      non-negative, and not all zero. The size of S is the size of P.
%
%    S = RESAMPSTR(P,M,N) returns an M by N matrix.
%
%    Default is to use no-sort resampling. For sorted resampling use
%      [PS,PI]=SORT(P);
%      S=PI(RESAMPSTR(PS));
%    Sorted re-sampling is slower but has slightly smaller
%    variance. Stratified resampling is unbiased, almost as fast as
%    deterministic resampling (RESAMPDET), and has only a slightly
%    larger variance.
%
%    In stratified resampling indices are sampled using random
%    numbers u_j~U[(j-1)/n,j/n], where n is length of P. Compare
%    this to simple random resampling where u_j~U[0,1].
%
%  Reference
%    Kitagawa, G., Monte Carlo Filter and Smoother for Non-Gaussian
%    Nonlinear State Space Models, Journal of Computational and
%    Graphical Statistics, 5(1):1-25, 1996.
%
%  See also 
%    RESAMPSIM, RESAMPRES, RESAMPDET

% Copyright (c) 2003-2004,2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin<2
  [m,n]=size(p);
elseif nargin==2
  n=m;
end
mn=m.*n;
pn=p./sum(p(:)).*mn;
s=zeros(m,n);
r=rand(1,mn);
k=0;
c=0;
for i=1:numel(pn)
  c=c+pn(i);
  if c>=1
    a=floor(c);
    c=c-a;
    s(k+[1:a])=i;
    k=k+a;
  end
  if k<mn && c>=r(k+1)
    c=c-1;
    k=k+1;
    s(k)=i;
  end
end
