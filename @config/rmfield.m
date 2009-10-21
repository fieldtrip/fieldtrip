function x = rmfield(x, key)

% RMFIELD Removes specified field from a CONFIGURATION object.

if isa(key, 'cell')
  for i=1:numel(key)
    x = rmfield(x, key{i});
  end
else
  x.value     = rmfield(x.value    , key);
  x.assign    = rmfield(x.assign   , key);
  x.reference = rmfield(x.reference, key);
  x.original  = rmfield(x.original , key);
end
