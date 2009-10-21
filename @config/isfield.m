function y = isfield(x, key)

% ISFIELD Get field names from a config object.

% all elements of an array would have the same fields, hence get the fields from the first element
y = isfield(x(1).value, key);

