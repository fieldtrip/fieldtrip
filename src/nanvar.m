function Y = nanvar(X, dim)

switch nargin
  case 1
    % MATLAB's var normalizes by n-1 when no dim is given. Adjust output
    % accordingly:
    Y = nanvar_base(X);
    n = nannumel(X);
    Y = Y .* (n ./ (n-1))
  case 2
    % In this case, the default of normalizing by n is expected.
    Y = nanvar_base(X, dim);    
  otherwise
    error ('Too many input arguments!')
end
