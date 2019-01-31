function Ainv = ft_inv(A, varargin)

% FT_INV computes a (regularized) inverse of a matrix.
%
% Use as
%  Ainv = ft_inv(A, ...)
% where optional additional arguments should be defined as key-value pairs.
%
% method    = 'tikhonov', 'svd', 'vanilla', 'moorepenrose', or 'winsorize',
%               the default method is 'tikhonov', and the matrix is regularized
%               with a weighted identity matrix before inversion, using the 
%               MATLAB pinv() function, using the lambda parameter.
%              'svd' results in a pseudoinverse
%               based on a truncated svd, fixing the number of singular values
%               according to tolerance/kappa (see below) before reassembling
%               the inverse. If the method specified is 'vanilla', a normal 
%               inv() is computed. If the method specified is 'moorepenrose',
%               a Moore-Penrose pseudoinverse is computed. If the method 
%               specified is 'winsorize', a truncated svd is computed, based
%               on the kappa/tolerance parameters, but in addition the 
%               singular values < lambda are replaced by the value according
%               to lambda.
% tolerance = scalar, reflects the fraction of the largest singular value
%               at which the singular value spectrum will be truncated. The
%               default value is 10*eps*max(size(A))
% kappa     = scalar integer, reflects the ordinal singular value at which
%               the singular value spectrum will be fixed. A value <=0 will
%               result in an all-zeros output matrix.
% lambda    = scalar, or string (expressed as a percentage), specifying the
%               regularization parameter for Tikhonov regularization, or 
%               the replacement value for winsorization. Lambda
%               specified as a percentage, e.g. '5%' will be converted into
%               a percentage of the average of trace(A).
% feedback  = boolean, false or true, to visualize the singular value spectrum
%               with the truncation level used.
% interactive = boolean, false or true, to manually specify a value for
%               kappa.
%
% kappa and tolerance are mutually exclusive, in case both are specified,
% tolerance takes precedence. Tolerance and kappa don't have an effect for
% the methods 'moorepenrose', 'tikhonov', 'vanilla' and 'winsorize'.  

% Copyright (C) 2019, DCCN, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.

method    = ft_getopt(varargin, 'method',    'tikhonov');
tolerance = ft_getopt(varargin, 'tolerance', []);
kappa     = ft_getopt(varargin, 'kappa',     []);
lambda    = ft_getopt(varargin, 'lambda',    0);
feedback  = istrue(ft_getopt(varargin, 'feedback', false));
interactive = istrue(ft_getopt(varargin, 'interactive', false));

[m, n] = size(A);
if m==0 || n==0
  % handle an empty input matrix
  Ainv = A';
  return;
end

if n < m
  % recurse into ft_inv on the transposed input
  Ainv = ft_inv(A',varargin{:})';
  return;
end

% perform an svd, because the handling of default tolerance etc may depend
% on this.
[U,S,V] = svd(A,0);
if m > 1
  s = diag(S);
elseif m == 1
  s = S(1);
else
  s = 0;
end

needkappa  = false;
needlambda = false;
switch method
  case 'svd'
    needkappa = true;
  case 'tikhonov'
    needlambda = true;
  case 'winsorize'
    needlambda = true;
    needkappa  = true;
end

if needlambda
  % translate the lambda argument
  if ischar(lambda) && lambda(end)=='%'
    ratio  = sscanf(lambda, '%f%%');
    ratio  = ratio/100;
    lambda = ratio * trace(A)./m;
  end
  
  if lambda~=0 && m~=n && strcmp(method, 'tikhonov')
    ft_error('with tikhonov regularization, the input matrix should be square');
  elseif lambda<0
    ft_error('a negative value for lambda is not allowed');
  end
end

if needkappa
  if isempty(kappa) && isempty(tolerance)
    % check whether the interacive mode is requested
    if interactive
      figure; plot(1:m, log10(s),'o-');
      kappa = input('Specify the number of dimensions at which the eigenvalue spectrum will be truncated: ');
      if isempty(kappa)
        kappa = m;
      end
    else
      % The default tolerance here is higher than the one used in MATLAB's
      % pinv, with a factor of 10. This was done by Robert Oostenveld in 2004,
      % to avoid numerical tolerance issues in the MATLAB version that was
      % used back then
      tolerance = 10 * max(m,n) * eps;
      kappa     = sum(s./s(1) > tolerance);
    end
  elseif ~isempty(kappa) &&  isempty(tolerance)
  elseif  isempty(kappa) && ~isempty(tolerance)
    kappa = sum(s./s(1) > tolerance);
  elseif ~isempty(kappa) && ~isempty(tolerance)
    ft_warning('both a kappa value and tolerance value are defined, using the kappa value for truncation');
  end
  
  % do a sanity check on the estimated kappa
  if kappa<=0 || kappa>m
    ft_warning('The kappa value falls outside the supported range, returning an all zeros matrix');
    Ainv = zeros(size(A'),class(A));
    return;
  end
end

switch method
  case 'vanilla'
    Ainv = inv(A);
  case 'svd'
    S    = diag(ones(kappa,1)./s(1:kappa));
    Ainv = V(:,1:kappa)*S*U(:,1:kappa)';
  case 'moorepenrose'
    Ainv = A'/(A*A');
  case 'tikhonov'
    % this is the way that fieldtrip handled covariance inversion.
    Ainv = pinv(A+eye(size(A,1)).*lambda, 10 * max(m,n) * s(1) * eps);
  case 'winsorize'
    % replace all singular values < lambda by lambda, and truncate
    s(s./s(1) < lambda) = lambda;
    S    = diag(ones(kappa,1)./s(1:kappa));
    Ainv = V(:,1:kappa)*S*U(:,1:kappa)';
  otherwise
    ft_error('Unsupported method for matrix inversion');
end

if feedback
  figure; hold on
  plot(1:m, log10(s),'o-');
  ylabel('singular values (log10)');
  abc = axis;
  if ~isempty(tolerance)
    plot([0 m], log10(tolerance).*[1 1], 'k');
  end
  if ~isempty(kappa)
    plot(kappa.*[1 1], abc([3 4]), 'r');
  end
end

