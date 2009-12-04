function y = fieldnames(x)

% FIELDNAMES Get field names from a config object.

% all elements of an array would have the same fields, hence get the fields from the first element
y = fieldnames(x(1).value);
