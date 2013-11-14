function [p,xx,sh]=kernelp(x,xx)
%KERNELP 1D Kernel density estimation of data, with automatic kernel width
%
%  [P,XX]=KERNELP(X,XX) return density estimates P in points XX,
%  given data and optionally ecvaluation points XX. Density
%  estimate is based on simple Gaussian kernel density estimate
%  where all kernels have equal width and this width is selected by
%  optimising plug-in partial predictive density. Works well with
%  reasonable sized X.
%

% Copyright (C) 2001-2003,2010 Aki Vehtari
%
% This software is distributed under the GNU General Public 
% Licence (version 3 or later); please refer to the file 
% Licence.txt, included with the software, for details.

if nargin < 1
  error('Too few arguments');
end
[n,m]=size(x);
if n>1 && m>1
  error('X must be a vector');
end
x=x(:);
[n,m]=size(x);

xa=min(x);xb=max(x);xab=xb-xa;
if nargin < 2
  nn=200;
  xa=xa-xab/20;xb=xb+xab/20;
  xx=linspace(xa,xb,nn);
else
  [mm,nn]=size(xx);
  if nn>1 && mm>1
    error('XX must be a vector');
  end
end
xx=xx(:);
m=length(x)/2;
rp=randperm(n);
xd=bsxfun(@minus,x(rp(1:m)),x(rp(m+1:end))');
sh=fminbnd(@(s) err(s,xd),xab/n*4,xab,optimset('TolX',xab/n*4));
p=mean(normpdf(bsxfun(@minus,x(rp(1:m)),xx'),0,sh));

function e=err(s,xd)
e=-sum(log(sum(normpdf(xd,0,s))));

function y = normpdf(x,mu,sigma)
y = -0.5 * ((x-mu)./sigma).^2 -log(sigma) -log(2*pi)/2;
y=exp(y);
