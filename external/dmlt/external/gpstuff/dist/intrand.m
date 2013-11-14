function r = intrand(p, m, n);
%INTRAND Random matrices from uniform integer distribution.
%
%   R = INTRAND(P) returns a matrix of random numbers chosen   
%   from the uniform integer distribution with parameter P. The 
%   minimum and maximum output integer numbers are P(1) and P(2). 
%   R = INTRAND(P,M,N) returns an M by N matrix.
%

% Copyright (c) 2000 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


if nargin < 1 | any(size([p(:)])~=[2 1])
  error('Incorrect range assignment in INTRAND');
end
if p(1)<p(2)
  rmin=ceil(p(1));
  rrange=floor(p(2))-rmin+1;
else
  rmin=ceil(p(2));
  rrange=floor(p(1))-rmin+1;
end

if nargin < 2
  r=floor(rand.*rrange+rmin);
else
  if nargin < 3
    n=m;
  end
  r=floor(rand([m n]).*rrange+rmin);
end
