function param = selparam(data)

% SELPARAM(DATA) extracts the fieldnames param of the structure data containing functional
% data, which have a dimensionality consistent with the dimord field in the data. Selparam
% is a helper function to selectdata

[dim, fn] = dimlength(data);
if numel(fn)==1 && strcmp(fn{1}, 'dimord'),
  %data is not source data
  dim = dim{1};
elseif numel(fn)>1,
  %data is source data new style
  ft_error('selparam only works with input data which is not of type ''source''');
end
ndim = length(dim);

% remove all 1's at the right side of dim
for i=ndim:-1:3
  if dim(i)==1,
    dim(i) = [];
    ndim   = ndim-1;
  else
    break
  end
end

fn = fieldnames(data);
sel = false(size(fn));
for i=1:numel(fn)
  siz    = size(data.(fn{i}));
  nsiz   = numel(siz);
  sel(i) = (nsiz==ndim) && all(siz==dim);
end
param = fn(sel);

% some fields should be excluded
param = setdiff(param, {'time', 'freq', 'channel','inside','label'});
