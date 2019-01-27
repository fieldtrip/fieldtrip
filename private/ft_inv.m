function Ainv = ft_inv(A, varargin)

% FT_INV computes a (regularized) inverse of a matrix.
%
% Use as
%  Ainv = ft_inv(A, ...)
% where optional additional arguments should be defined as key-value pairs.
%
% method    = 'svd', 'vanilla', 'moorepenrose' or 'winsorize', the default
%               method is svd, which results in a pseudoinverse based on
%               an svd, fixing the number of singular values according to
%               tolerance/kappa (see below) before reassembling the inverse.
%               If the method specified is 'vanilla', a normal inv() is
%               computed. If the method specified is 'moorepenrose' is
%               specified, a Moore-Penrose pseudoinverse is computed. If
%               the method specified is 'winsorize' the singular values s are
%               clipped according to tolerance/kappa, but the original
%               number of singular values is used for reassembling the
%               inverse.
% tolerance = scalar, reflects the fraction of the largest singular value
%               at which the singular value spectrum will be truncated. The
%               default value is 10*eps*max(size(A))
% kappa     = scalar integer, reflects the ordinal singular value at which
%               the singular value spectrum will be fixed. A value <=0 will
%               result in an all-zeros output matrix.
% lambda    = scalar, or string (expressed as a percentage), specifying the
%               regularization parameter for diagonal loading. Lambda
%               specified as a percentage will be converted into a
%               percentage of the average of trace(A).
% feedback  = boolean, false or true, to visualize the singular value spectrum
%               with the truncation level used.
% interactive = boolean, false or true, to manually specify a value for
%               kappa.
%
% kappa and tolerance are mutually exclusive, in case both are specified,
% tolerance takes precedence. Tolerance and kappa don't have an effect for
% the methods 'moorepenrose' and 'vanilla'

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

method    = ft_getopt(varargin, 'method', 'svd');
tolerance = ft_getopt(varargin, 'tolerance', []);
kappa     = ft_getopt(varargin, 'kappa', []);
lambda    = ft_getopt(varargin, 'lambda', 0);
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

% translate the lambda argument
if ischar(lambda) && lambda(end)=='%'
  ratio  = sscanf(lambda, '%f%%');
  ratio  = ratio/100;
  lambda = ratio * trace(A)./m;
end

if lambda~=0 && m~=n
  ft_error('with lambda regularization, the input matrix should be square');
end

if m==n && lambda>0
  % regularize
  A = A + eye(m).*lambda;
elseif m~=n && lambda~=0
  ft_error('when lambda larger than zero, the input matrix must be square');
elseif m==n && lambda<0
  ft_error('a negative value for lambda is not allowed');
end

% perform an svd
[U,S,V] = svd(A,0);
if m > 1
  s = diag(S);
elseif m == 1
  s = S(1);
else
  s = 0;
end

if isempty(kappa) && isempty(tolerance)
  % check whether the interacive mode is requested
  if interactive
    figure; plot(1:m, log10(s),'o-');
    kappa = input('Specify the number of dimensions at which the eigenvalue spectrum will be truncated: ');
    if isempty(kappa)
      kappa = m;
    end
    usekappa     = true;
    usetolerance = false;
  else
    % The default tolerance here is higher than the one used in MATLAB's
    % pinv, with a factor of 10. This was done by Robert Oostenveld in 2004,
    % to avoid numerical tolerance issues in the MATLAB version that was
    % used back then
    tolerance = 10 * max(m,n) * eps;
    kappa     = sum(s./s(1) > tolerance);
    
    usekappa     = false;
    usetolerance = true;
  end
elseif ~isempty(kappa) &&  isempty(tolerance)
  usekappa     = true;
  usetolerance = false;
elseif  isempty(kappa) && ~isempty(tolerance)
  usekappa     = false;
  usetolerance = true; % needed for winsorization
  kappa = sum(s./s(1) > tolerance);
elseif ~isempty(kappa) && ~isempty(tolerance)
  ft_warning('both a kappa value and tolerance value are defined, using the kappa value for truncation');
  usekappa     = true;
  usetolerance = false;
end

% do a sanity check on the estimated kappa
if kappa<=0 || kappa>m
  ft_warning('The kappa value falls outside the supported range, returning an all zeros matrix');
  Ainv = zeros(size(A'),class(A));
  return;
end

switch method
  case 'vanilla'
    Ainv = inv(A);
  case 'svd'
    S    = diag(ones(kappa,1)./s(1:kappa));
    Ainv = V(:,1:kappa)*S*U(:,1:kappa)';
  case 'moorepenrose'
    Ainv = A'/(A*A');
  case 'winsorize'
    % replace all singular values < the kappa'th singular by the kappa'th
    % singular value
    if usekappa
      s = max(s,s(kappa));
    elseif usetolerance
      s(s./s(1) < tolerance) = tolerance;
    end
    S    = diag(1./s);
    Ainv = V*S*U';
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
  plot(kappa.*[1 1], abc([3 4]), 'r');
end

