% NANVAR provides a replacement for MATLAB's nanvar that is almost
% compatible.
%
% For usage see VAR. Note that the weight-vector is not supported. If you
% need it, please file a ticket at our bugtracker.

function Y = nanvar(X, dim)

switch nargin
  case 1
    % MATLAB's var normalizes by n-1 when no dim is given. Adjust output
    % accordingly:
    Y = nanvar_base(X);
    n = nannumel(X);
    Y = Y ./ (n-1);
  case 2
    % In this case, the default of normalizing by n is expected.
    Y = nanvar_base(X, dim);    
    n = nannumel(X);
    Y = Y ./ n;
  otherwise
    error ('Too many input arguments!')
end
