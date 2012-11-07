% NANVAR provides a replacement for MATLAB's nanvar that is almost
% compatible.
%
% For usage see VAR. Note that the weight-vector is not supported. If you
% need it, please file a ticket at our bugtracker.

function Y = nanvar(X, w, dim)

switch nargin
  case 1
    % VAR(x)
    % Normalize by n-1 when no dim is given. 
    Y = nanvar_base(X);
    n = nannumel(X);
    w = 0;

  case 2
    % VAR(x, 1)
    % VAR(x, w)
    % In this case, the default of normalizing by n is expected.
    Y = nanvar_base(X);
    n = nannumel(X);

  case 3
    % VAR(x, w, dim)
    % if w=0 normalize by n-1, if w=1 normalize by n.
    Y = nanvar_base(X, dim);    
    n = nannumel(X, dim);

  otherwise
    error ('Too many input arguments!')
end

% Handle different forms of normalization:
if numel(w) == 0
  % empty weights vector defaults to 0
  w = 0;
end

if numel(w) ~= 1
  error('Weighting vector w is not implemented! Please file a bug.');
end

if ~isreal(X)
  Y = real(Y) + imag(Y);
  n = real(n);
end

if w == 1
  Y = Y ./ n;
end

if w == 0
  Y = Y ./ max(1, (n - 1));  % don't divide by zero!
end
