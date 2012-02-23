function [volume, permutevec] = volumepermute(volume, permutevec)

% do a low-level check on the input data
if ~isfield(volume, 'transform'), error('the input volume needs a transformation matrix'); end
if ~isfield(volume, 'dim'),       error('the input volume needs a dim field');             end

if isfield(volume, 'outside')
end

% define some variable locally
T   = volume.transform;
dim = volume.dim;

% determine which fields can be permuted
fnames = fieldnames(volume);
sel    = false(1,numel(fnames));
for k = 1:numel(fnames)
  tmp = volume.(fnames{k});
  if isnumeric(tmp) && numel(tmp)==prod(dim) && ndims(tmp)==3
    sel(k) = true;
  end
end
fnames = fnames(sel);

if nargin<2
  permutevec = 'auto';
end

if ischar(permutevec)
  % determine the order of permutation to make the transformatiom matrix approximately diagonal
  [dum, m1]  = max(abs(T(1,1:3)));
  [dum, m2]  = max(abs(T(2,1:3)));
  [dum, m3]  = max(abs(T(3,1:3)));
  permutevec = [m1 m2 m3]; 
end

if ~all(permutevec==[1 2 3])
  % do the permutation on the numeric data
  for k = 1:numel(fnames)
    volume = setsubfield(volume, fnames{k}, permute(getsubfield(volume, fnames{k}), permutevec));
  end
  
  % do the permutation on the transformation matrix and dim
  volume.transform = volume.transform(:,[permutevec 4]);
  volume.dim       = volume.dim(permutevec);
else
  % do nothing
end