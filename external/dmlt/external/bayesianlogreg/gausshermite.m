function [x,w] = gausshermite(n);

% In:
%
% n: number of abscissas
%
% Out:
%
% x,w: abscissas and weights
%
% Taken from numerical recipes but now such that
%
% int dx exp(-x^2/2) / sqrt(2 pi) f(x) = sum_i w_i f(x_i)

if nargin < 1, n = []; end
if ~length(n), n=9; end

m = (n+1)/2;
x = zeros(n,1);
w = zeros(n,1);

for i=1:m,

  if i == 1,
    z = sqrt(2*n+1)-1.85575*(2*n+1)^(-0.16667);
  elseif i == 2,
    z = z - 1.14*n^0.426/z;
  elseif i == 3,
    z = 1.86*z-0.86*x(1);
  elseif i == 4,
    z = 1.91*z-0.91*x(2);
  else
    z = 2.0*z-x(i-2);
  end
  for its=1:10,
    p1 = 1/sqrt(sqrt(pi));
    p2 = 0.0;
    for j=1:n,
      p3 = p2;
      p2 = p1;
      p1 = z*sqrt(2/j)*p2-sqrt((j-1)/j)*p3;
    end
    pp = sqrt(2*n)*p2;
    z1 = z;
    z = z1-p1/pp;
    if (abs(z-z1) < eps), break, end;
  end
  x(i) = z;
  x(n+1-i) = -z;
  w(i) = 2/(pp*pp);
  w(n+1-i) = w(i);

end

% this x and w for
% int dx exp(-x^2) f(x) = sum_i w_i f(x_i)
% now transform to
% int dx exp(-x^2/2) / sqrt(2 pi) f(x) = sum_i w_i f(x_i)

w = w/sqrt(pi);
x = x*sqrt(2);
