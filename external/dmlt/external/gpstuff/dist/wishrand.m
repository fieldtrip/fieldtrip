function R = wishrand(S,nu);
%WISHRND Random matrices from Wishard distribution.
%   R = WISHRAND(S,N) returns a matrix of random numbers chosen   
%   from the central Wishard distribution with parameters S and NU.
%
%   S is a symmetric positive definite scale matrix
%   NU is degrees of freedom
%
%   Note: E[R]=S
%
%   References: Gelman, Carlin, Stern, & Rubin (1995), Bayesian Data Analysis,
%                 Chapman & Hall, pp. 474, 478, 480-481.
%               Gentle (2003), Random Number Generation and Monte Carlo
%                 Methods, 2nd ed, Springer, p. 199, Algorithm 5.8. 
%               Smith & Hocking (1972), Algorithm AS 53: Wishard Variate
%                 Generator, Applied Statistics, 21(3), pp. 341-345.
%
%	See also INVWISHRAND

% Copyright (c) 1999-2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 2
  error('Requires two input arguments.'); 
end;

[d d2] = size(S);
if d ~= d2
  error('Matrix S must be square');
end

% Algorithm 5.8 step 0.
[T,p]=chol(S);
if p > 0
  error('Matrix S must be positive definite.');
end

if (nu >= d) && (nu == round(nu))
  % distribution is proper and degrees of freedom is integer
  % brute-force may be used, which is surprisingly faster in Matlab
  % at least up to d>10, nu>1000
  % Algorithm described e.g. in (Gelman et al., 1995)
  Y = T'*randn(d,nu);
  R = Y*Y'./nu;
else
  % distribution is not proper or degrees of freedom is not integer
  % (Gentle, 2003) Algorithm 5.8
  % Algorithm 5.8 step 1.
  Z=zeros(d);
  for j=2:d
    Z(1:j,j)=randn(j,1);
  end
  % Algorithm 5.8 step 2. Note that there is error in the book.
  % Book says nu-i, while it should be nu+1-i
  % See errata <http://www.scs.gmu.edu/~jgentle/rngbk/errata.htm>
  % gamrand below is same as chi2rnd(nu-(1:d)+1), but much faster using c-code
  y=gamrand((nu+1-(1:d)),nu+1-(1:d)); 
  % Algorithm 5.8 step 3.
  Z2=Z.^2;
  sy=sqrt(y);
  B=zeros(d);
  B(1,1)=y(1);
  % In Matlab 6.5 following loop is accelerated and thus quite fast
  for j=2:d
    B(j,j)=y(j)+sum(Z2(1:(j-1),j));
    B(1,j)=Z(1,j).*sy(1);
    B(j,1)=B(1,j);
    for i=2:(j-1)
      B(i,j)=Z(i,j).*sy(i)+sum(Z(1:(i-1),i).*Z(1:(i-1),j));
      B(j,i)=B(i,j);
    end
  end
  % Algorithm 5.8 step 4.
  R=T'*B*T/nu;
end
