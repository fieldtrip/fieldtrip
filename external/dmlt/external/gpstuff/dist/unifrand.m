function r = unifrand(a,b,m,n)
% UNIFRAND - Generate unifrom random numberm from interval [A,B]
%
%    R = UNIFRAND(MU) returns a matrix of random numbers chosen   
%    from the uniform distribution with parameters A and B.
%    The size of R is the common size of A and B.
%    Alternatively, R = UNIFRAND(A,B,M,N) returns an M by N matrix. 
% 

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if nargin<2
  error('Not enough arguments')
end
if a>b
  error('A must be smaller than B')
end
if nargin<3
  if not(all(size(a)==size(b)))
    error('A and B must be same size');
  end
  [m,n]=size(a);
end
if nargin<4
  n=m;
end
d=b-a;
r=rand(m,n)*d+a;
