function s = resampdet(p,m,n);
%RESAMPDET Deterministic resampling
%
%   Description
%   S = RESAMPDET(P) returns a new set of indices according to the
%   probabilities P. P is array of probabilities, which are not
%   necessarily normalized, though they must be non-negative, and
%   not all zero. The size of S is the size of P. 
%
%   S = RESAMPDET(P,M,N) returns M by N matrix.
%
%   S = RESAMPDET(P,M) returns M by M matrix.
%
%   Default is to use no-sort resampling. For sorted resampling use
%    [PS,PI]=SORT(P);
%    S=PI(RESAMPDET(PS));
%   Sorted re-sampling is slower but has smaller variance. Note
%   that deterministic resampling is not unbiased. Stratified
%   resampling (RESAMPSTR) is unbiased, almost as fast as
%   deterministic resampling, and has only slightly larger
%   variance.
%
%   In deterministic resampling indices are sampled using
%   deterministic numbers u_j~(j-a)/n, for fixed a in [0,1) and
%   n is length of P. Compare this to simple random resampling
%   where u_j~U[0,n]. See, Kitagawa, G., Monte Carlo Filter and
%   Smoother for Non-Gaussian Nonlinear State Space Models,
%   Journal of Computational and Graphical Statistics, 5(1):1-25,
%   1996. 
%
%   See also RESAMPSIM, RESAMPRES, RESAMPSTR

% Copyright (c) 2003-2004 Aki Vehtari

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
fpn=floor(pn);
s=zeros(m,n);
k=0;
c=0.5;
for i=1:numel(p)
  if pn(i)>=1
    a=fpn(i);
    pn(i)=pn(i)-a;
    s(k+[1:a])=i;
    k=k+a;
  end
  c=c+pn(i);
  if c>=1
    k=k+1;
    s(k)=i;
    c=c-1;
  end
end
