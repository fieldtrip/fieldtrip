function x = set(x, key, val)

% SET Assign a new value to the field of a config object.

if ~isfield(x.value, key)
  % initialize the counters for this field
  % see the explaination about side effects of the increment function in config.m
  x.assign.(key)    = deepcopy(0); % ensure that a unique scalar is created for each counter
  x.reference.(key) = deepcopy(0); % ensure that a unique scalar is created for each counter
  x.original.(key)  = deepcopy(0); % ensure that a unique scalar is created for each counter
end

x.value.(key) = val;
increment(x.assign.(key));
