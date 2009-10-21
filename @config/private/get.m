function y = get(x, key, inc)

% GET Return the value of a field in a config object.

if nargin<3
  % the default is to increment the reference counter, an exception is when
  % a structure has to be recursed to assign a value to a substructure
  inc = true;
end

if isfield(x.value, key)
  y = x.value.(key);
  if inc
    increment(x.reference.(key));
  end
else
  error(sprintf('Reference to non-existent field ''%s''.', key));
end
