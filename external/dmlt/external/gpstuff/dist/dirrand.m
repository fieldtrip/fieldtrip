function r = dirrand(m,n)
% DIRRAND - Uniform Dirichlet random vectors
%
%    R = DIRRAND(M) returns vector of length M chosen from
%    Dirichlet(1,...,1) distribution.
%    R = DIRRAND(M,N) returns N such vectors in MxN matrix.
%
%   References: Gelman, Carlin, Stern, & Rubin (2003), Bayesian Data Analysis,
%               Chapman & Hall, pp. 582.

% Copyright (c) 1998-2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin < 2
  n=1;
end

% Generate random numbers from Gamma(1,1) distribution...
r=rand(m,n);
r=-reallog(r);
% ...and scale sum to unity
rs=sum(r);
for i1=1:m
  r(i1,:)=r(i1,:)./rs;
end
