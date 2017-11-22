function result = mergecellstruct(result)

% MERGECELLSTRUCT is a helper function for FT_TEST

if iscell(result)
  if numel(result)>1
    result = appendstruct(result{:});
  else
    result = result{1};
  end
end
