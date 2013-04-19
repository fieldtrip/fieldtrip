function lp=dir_lpdf(x,a)
%DIR_LPDF   Log probability density function of uniform Dirichlet
%           distribution
%
%       Description:
%       LP = DIR_LPDF(X, A) returns the lpdf of Dirichlet distribution 
%       with A at X

% Copyright (c) 2000 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


lp=gammaln(sum(a))-sum(gammaln)+sum(log(x).^(a-1));
