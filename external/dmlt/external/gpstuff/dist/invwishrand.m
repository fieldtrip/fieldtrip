function r = invwishrand(S,nu);
%INVWISHRND Random matrices from inverse Wishard distribution.
%
%   R = INVWISHRAND(S,N) returns a matrix of random numbers chosen   
%   from the central inverse Wishard distribution with parameters S and NU.
%
%   S is a symmetric positive definite scale matrix
%   NU is degrees of freedom
%
%   Note: E[R]=S*nu/(nu-k-1)
%
%	See also WISHRAND

% Copyright (c) 1999 Aki Vehtari

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

[t,p]=chol(S);
if p < 0
  error('Matrix S must be positive definite.');
end

r=inv(wishrand(inv(S),nu));
