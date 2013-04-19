function brnd = betarand(a,b,r,c)
%  BETARAND Random matrices from beta distribution.
%     R = BETARAND(A,B) returns an array of random numbers chosen from the
%     beta distribution with parameters A and B.  The size of R is the common
%     size of A and B if both are arrays.  If either parameter is a scalar,
%     the size of R is the size of the other parameter.
%
%     R = BETARAND(A,B,M,N) returns an M-by-N matrix.

% Copyright (c) 1995-1997, 2005-2007 Kurt Hornik

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
if (nargin > 1)
  if (~isscalar(a) || ~isscalar(b))
    if (size(a,1) ~= size(b,1))
      error ('betarnd: a and b must be of common size or scalar');
    end
  end
end

if (nargin == 4)
  if (~(isscalar(r) && (r > 0) && (r == round(r))))
    error ('betarnd: r must be a positive integer');
  end
  if (~(isscalar (c) && (c > 0) && (c == round (c))))
    error ('betarnd: c must be a positive integer');
  end
  sz = [r, c];
  
  if (any (size (a) ~= 1) && (length (size (a)) ~= length (sz) || any (size (a) ~= sz)))
    error ('betarnd: a and b must be scalar or of size [r,c]');
  end
elseif (nargin == 3)
  if (isscalar (r) && (r > 0))
    sz = [r, r];
  elseif (isvector(r) && all (r > 0))
    sz = r(:)';
  else
    error ('betarnd: r must be a positive integer or vector');
  end
  
  if (any (size (a) ~= 1) && (length (size (a)) ~= length (sz) || any (size (a) ~= sz)))
    error ('betarnd: a and b must be scalar or of size sz');
  end
elseif (nargin == 2)
  sz = size(a);
else
  error('must provide atleast 2 parameters')
end

if (isscalar(a) && isscalar(b))
  if (find ((a < 0) | isinf(a) | (b < 0) | isinf(b)))
    brnd = NaN * ones (sz);
  else
    r1 = gamrand(a,2.*a,sz(1),sz(2));
    brnd = r1 ./ (r1 + gamrand(b,2.*b,sz(1),sz(2)));
  end
else
  brnd = zeros (sz);
  
  k = find ((a < 0) | isinf(a) | (b < 0) | isinf(b));
  if (any (k))
    brnd(k) = NaN * ones (size (k));
  end
  
  k = find ((a > 0) & (a < Inf) & (b > 0) & (b < Inf));
  if (any (k))
    r1 = gamrnd(a(k),1,size(k,1),size(k,2));
    brnd(k) = r1 ./ (r1 + gamrnd(b(k),1,size(k,1),size(k,2)));
  end
end

end