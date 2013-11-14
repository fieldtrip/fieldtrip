function r = normtzrand(mu,s,lr)
%NORMZTRAND Random draws from a normal distribution truncated by zero
%
%    R = NORMTZRAND(MU,S,LR) returns a matrix of random numbers
%    chosen from the normal distribution truncated by zero with
%    mean MU, std S, and whether truncation is to left or right by
%    LR, where left=1 and right=0. The size of R is the common size
%    of MU, S and LR.

% Copyright (c) 2003 Aki Vehtari

% Implementation based on code by Jim Albert, see
% Ordinal Data Modeling by Valen Johnson and James Albert
% Springer-Verlag, New York, 1999.

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin ~= 3
  error('Requires four input arguments.'); 
end;

lr=logical(lr);
a = 0.5*(1+erf(-mu./s./sqrt(2)));
a(lr) = 1-a(lr);
u = rand(size(mu)).*a;
r = sqrt(2).*erfinv(2*u-1).*s;
r(lr)=-r(lr);
r = r + mu;
