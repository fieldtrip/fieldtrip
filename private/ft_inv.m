function Y = ft_inv(X, varargin)

% FT_INV computes a matrix inverse with optional regularization.
%
% Use as
%  Y = ft_inv(X, ...)
%
% Additional options should be specified in key-value pairs and can be
%   method    = string, method for inversion and regularization (see below).
%               The default method is 'lavrentiev'.
%   lambda    = scalar value, or string (expressed as a percentage), specifying 
%               the regularization parameter for Lavrentiev or Tikhonov 
%               regularization, or the replacement value for winsorization. 
%               When lambda is specified as a string containing a percentage, 
%               e.g. '5%', it will be computed as the percentage of the average 
%               eigenvalue.
%   kappa     = scalar integer, reflects the ordinal singular value at which
%               the singular value spectrum will be truncated.
%   tolerance = scalar, reflects the fraction of the largest singular value
%               at which the singular value spectrum will be truncated.
%               The default is 10*eps*max(size(X)).
%   feedback  = boolean, to visualize the singular value spectrum with the 
%               lambda regularization and kappa truncation.
%
% The supported methods are:
%
% 'vanilla' - the MATLAB inv() function is used for inversion, no regularization is
% applied.
%
% 'moorepenrose' - the Moore-Penrose pseudoinverse is computed, no regularization is
% applied.
%
% 'tsvd' - this results in a pseudoinverse based on a singular value decomposition,
% truncating the singular values according to either kappa or tolerance parameter
% before reassembling the inverse.
%
% 'tikhonov' - the matrix is regularized according to the Tikhonov method using the
% labmda parameter, after which the truncated svd method (i.e. similar to MATLAB
% pinv) is used for inversion.
%
% 'lavrentiev' - the matrix is regularized according to the Lavrentiev method with a
% weighted identity matrix using the labmda parameter, after which the truncated svd
% method (i.e. similar to MATLAB pinv) is used for inversion.
%
% 'winsorize' - a truncated svd is computed, based on either kappa or tolerance
% parameters, but in addition the singular values smaller than lambda are replaced by
% the value according to lambda.
%
% Both for the lambda and the kappa option you can specify 'interactive' to pop up an
% interactive display of the singular value spectrum that allows you to click in the figure. 
%
% Rather than specifying kappa, you can also specify the tolerance as the ratio of
% the largest eigenvalue at which eigenvalues will be truncated.
%
% See also INV, PINV, CONDEST, RANK

% Copyright (C) 2019, DCCN, Jan-Mathijs Schoffelen, Robert Oostenveld
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

method      = ft_getopt(varargin, 'method',    'lavrentiev');
tolerance   = ft_getopt(varargin, 'tolerance', []);
kappa       = ft_getopt(varargin, 'kappa',     []);
lambda      = ft_getopt(varargin, 'lambda',    0);
feedback    = istrue(ft_getopt(varargin, 'feedback', false));

[m, n] = size(X);
if m==0 || n==0
  % handle an empty input matrix
  Y = X';
  return;
end

switch method
  case 'vanilla'
    needkappa  = false;
    needlambda = false;
  case 'moorepenrose'
    needkappa  = false;
    needlambda = false;
  case 'tsvd'
    needkappa  = true;
    needlambda = false;
  case 'tikhonov'
    needkappa  = true;
    needlambda = true;
  case 'lavrentiev'
    needkappa  = true;
    needlambda = true;
  case 'winsorize'
    needlambda = true;
    needkappa  = true;
  otherwise
    ft_error('unsupported method "%s"', method);
end

if needlambda || needkappa || feedback
  % perform an svd, this is used to determine the default parameters and for most inverse methods
  [U,S,V] = svd(X,0);
  if m > 1
    s = diag(S);
  elseif m == 1
    s = S(1);
  else
    s = 0;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambda is the regularization factor that inflates the eigenvalues

if needlambda
  if ischar(lambda) && lambda(end)=='%'
    % translate the lambda argument into a scalar
    ratio  = sscanf(lambda, '%f%%');
    ratio  = ratio/100;
    lambda = ratio * mean(s);
  elseif ischar(lambda) && strcmp(lambda, 'interactive')
    figure
    semilogy(1:m, s,'o-');
    title('Please specify lambda by clicking in the figure');
    fprintf('Please specify lambda by clicking in the figure\n');
    [x,lambda] = ginput(1);
    % compute the difference between the vertical position that was clicked and the corresponding eigenvalue
    lambda = lambda - s(round(x));
    if lambda<0
      lambda = 0;
    end
    if isempty(lambda)
      lambda = 0;
    end
  end
  
  assert(lambda>=0, 'a negative value for lambda is not allowed');
else
  assert(isempty(lambda) || lambda==0, 'this method is incompatible with the specification of lambda');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kappa is the truncation parameter, expressed as an integer
% Alternatively the tolerance can be specified, expressed in matrix units

if needkappa
  if ischar(kappa) && strcmp(kappa, 'interactive')
    figure
    semilogy(1:m, s,'o-');
    title('Please specify kappa by clicking in the figure');
    fprintf('Please specify kappa by clicking in the figure.');
    [kappa,~] = ginput(1);
    if kappa<0
      kappa = 0;
    end
    if isempty(kappa) || kappa>m
      kappa = m;
    end
  end
  
  if isempty(kappa) &&  isempty(tolerance)
    % The default tolerance here is higher than the one used in MATLAB's pinv, with a
    % factor of 10. This was done by Robert Oostenveld in 2004, to avoid numerical
    % tolerance issues in the MATLAB version that was used back then.
    tolerance = 10 * max(m,n) * eps;
    kappa     = sum(s./s(1) > tolerance);
  elseif  isempty(kappa) && ~isempty(tolerance)
    kappa = sum(s./s(1) > tolerance);
  elseif ~isempty(kappa) &&  isempty(tolerance)
    assert(kappa>=0 && kappa<=m, 'invalid specification of kappa')
  elseif ~isempty(kappa) && ~isempty(tolerance)
    ft_warning('both kappa and tolerance value are defined, using the kappa value for truncation');
  end
  
  assert(kappa>=0, 'a negative value for kappa is not allowed');
else
  assert(isempty(kappa) || kappa==m, 'this method is incompatible with the specification of kappa');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% provide initial feedback

if feedback
  % plot the initial singular values
  figure
  semilogy(1:m, s, 'o-');
  hold on
  ylabel('singular values');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the actual inverse

switch method
  case 'vanilla'
    Y = inv(X);
  case 'moorepenrose'
    Y = X'/(X*X');
  case 'tsvd'
    S = diag(1./s(1:kappa));
    Y = V(:,1:kappa)*S*U(:,1:kappa)';
  case 'tikhonov'
    % truncate singular values at kappa, inflate all non-zero singular values with lambda
    % see https://www.uni-muenster.de/AMM/num/Vorlesungen/IP_WS07/skript.pdf example 3.15
    S = diag(s(1:kappa) ./ (s(1:kappa).^2 + lambda));
    if feedback
      semilogy(1:kappa, 1./diag(S), 'k+');
    end
    Y = V(:,1:kappa)*S*U(:,1:kappa)';
  case 'lavrentiev'
    % truncate singular values at kappa, inflate all non-zero singular values with lambda
    % see https://www.uni-muenster.de/AMM/num/Vorlesungen/IP_WS07/skript.pdf example 3.14
    S = diag(1 ./ (s(1:kappa) + lambda));
    if feedback
      semilogy(1:kappa, 1./diag(S), 'k+');
    end
    Y = V(:,1:kappa)*S*U(:,1:kappa)';
  case 'winsorize'
    % replace all singular values smaller than lambda by lambda, and truncate by kappa
    s(s<lambda) = lambda;
    S = diag(ones(kappa,1)./s(1:kappa));
    if feedback
      semilogy(1:kappa, 1./diag(S), 'k+');
    end
    Y = V(:,1:kappa)*S*U(:,1:kappa)';
  otherwise
    ft_error('unsupported method "%s"', method);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize the feedback

if feedback
  abc = axis;
  if ~isempty(tolerance) && tolerance>0
    % plot a horizontal line indicating the tolerance
    plot(abc([1 2]), tolerance.*[1 1], 'g');
  end
  if ~isempty(kappa) && kappa<m
    % plot a vertical line indicating kappa truncation
    plot(kappa.*[1 1], abc([3 4]), 'r');
  end
  axis(abc);
end
