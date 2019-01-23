function Ainv = ft_inv(A, varargin)

% FT_INV computes a regularized inverse of a matrix.
%
% Use as 
%  Ainv = ft_inv(A, ...)
% where optional additional arguments should be defined as key-value pairs.
%
% method    = 'pinv', 'inv', or 'winsorize'
% tolerance = scalar, blabla
% kappa     = scalar integer
% lambda    = scalar, or string 
% feedback  = boolean
%
% kappa and tolerance are mutually exclusive, in case both are specified,
% tolerance takes precedence

method    = ft_getopt(varargin, 'method', 'pinv');
tolerance = ft_getopt(varargin, 'tolerance', []);
kappa     = ft_getopt(varargin, 'kappa', []);
lambda    = ft_getopt(varargin, 'lambda', 0);
feedback  = istrue(ft_getopt(varargin, 'feedback', false));
interactive = istrue(ft_getopt(varargin, 'interactive', false));

[m, n] = size(A);

if n < m
  % recurse into ft_inv on the transposed input
  Ainv = ft_inv(A',varargin{:})';
  return
elseif m <= n 
  % here the work starts
  if ~isempty(kappa) && kappa > m
    ft_error('the kappa value for truncation is larger than the dimensionality of the matrix');
  end
  
  if ~isempty(kappa) && ~isempty(tolerance)
    ft_warning('both a kappa value and tolerance value are defined, using the kappa value for truncation');
  end
  
  if isempty(lambda)
    lambda = 0;
  elseif ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
    ratio  = sscanf(lambda, '%f%%');
    ratio  = ratio/100;
    lambda = ratio * trace(A)./m;
  end
  
  if lambda~=0 && m~=n
    ft_error('with lambda regularization, the input matrix should be square');
  end
  
  if m==n && lambda~=0
    % regularize
    A = A + eye(m).*lambda;
  elseif lambda~=0
    ft_error('regularization is not possible when the matrix is not square');
  end
  
  [U,S,V] = svd(A,0);
  if m > 1
    s = diag(S);
  elseif m == 1
    s = S(1);
  else
    s = 0;
  end
  
  if isempty(kappa) && isempty(tolerance) && interactive
    figure; plot(1:m, log10(s),'o-');
    kappa = input('Specify the number of dimensions at which the eigenvalue spectrum will be truncated: ');
    if isempty(kappa)
      kappa = m;
    end
  end
    
  elseif isempty(kappa) && isempty(tolerance)
    % put the tolerance for the pinv to its default 
    tolerance = 10 * max(m,n) * max(s) * eps;
  end
  
  if ~isempty(kappa)
    % use kappa
  elseif ~isempty(tolerance)
    kappa = sum(s > tolerance);
  end
  
  if feedback
    figure; hold on 
    plot(1:m, log10(s),'o-');
    ylabel('singular values (log10)');
    abc = axis;
    if ~isempty(tolerance)
      plot([0 m], log10(tolerance).*[1 1], 'k');
    end
    plot(kappa.*[1 1], abc([3 4]), 'r');
  end
  
  if kappa==0
    Ainv = zeros(size(A'),class(A));
    return;
  end
  
  switch method
    case 'inv'
      Ainv = inv(A);
    case 'pinv'
      s    = diag(ones(kappa,1)./s(1:kappa));
      Ainv = V(:,1:kappa)*s*U(:,1:kappa)';
    case 'winsorize'
    otherwise 
      ft_error('Unsupported method for matrix inversion');
  end
  
end
