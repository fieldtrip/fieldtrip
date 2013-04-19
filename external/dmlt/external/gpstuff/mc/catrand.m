function r = catrand(p, m, n);
%CATRAND Random matrices from categorical distribution.
%
%  Description
%    R = CATRAND(P) returns a matrix of random numbers chosen from
%      the categorical distribution with parameter P. P is array of
%      probabilities, which are not necessarily normalized, though
%      they must be non-negative, and not all zero. The size of R
%      is the size of P.
%
%    R = CATRAND(P,M,N) returns an M by N matrix.

% Copyright (c) 1999-2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


p=p(:);
pc=cumsum(p);
pc=pc./pc(end);
if nargin < 2
  m=1;n=1;
elseif nargin <3
  n=m;
end
r=binsgeq(pc,rand(m,n));
