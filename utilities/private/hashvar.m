function hash = hashvar(var)

% HASHVAR computes a hash on the input variable, after removing any
% FieldTrip specific fields that are not considered to be part of the
% data.

if isstruct(var) && isfield(var, 'cfg')
  var = rmfield(var, 'cfg');
end

hash = CalcMD5(mxSerialize(var));

