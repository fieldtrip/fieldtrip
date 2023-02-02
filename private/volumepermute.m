function [volume, permutevec] = volumepermute(volume, permutevec)

% VOLUMEPERMUTE
%
% See also VOLUMEFLIP, ALIGN_IJK2XYZ, ALIGN_XYZ2IJK

% do a low-level check on the input data
if ~isfield(volume, 'transform'), ft_error('the input volume needs a transformation matrix'); end
if ~isfield(volume, 'dim'),       ft_error('the input volume needs a dim field');             end

% define some variable locally
T   = volume.transform;
dim = volume.dim;

% determine which fields can be permuted
fn = parameterselection('all', volume);

if nargin<2
  permutevec = 'auto';
end

if isequal(permutevec, 'auto')
  % determine the order of permutation to make the transformatiom matrix approximately diagonal
  [dum, m1]  = max(abs(T(1,1:3)));
  [dum, m2]  = max(abs(T(2,1:3)));
  % [dum, m3]  = max(abs(T(3,1:3)));
  [dum, m3]  = setdiff(1:3, [m1 m2]); % whichever dimension remains
  permutevec = [m1 m2 m3];
end

if ~all(permutevec==[1 2 3])
  % do the permutation on the numeric data
  for k = 1:numel(fn)
    tmp = volume.(fn{k});
    switch ndims(tmp)
      case 3
        volume.(fn{k}) = permute(tmp, permutevec);
      case 4
        volume.(fn{k}) = permute(tmp, [permutevec 4]);
      case 5
        volume.(fn{k}) = permute(tmp, [permutevec 4 5]);
      otherwise
        ft_error('cannot permute %s', fn{k});
    end
  end

  % do the permutation on the transformation matrix and dim
  volume.transform = volume.transform(:,[permutevec 4]);
  volume.dim       = volume.dim(permutevec);
else
  % do nothing
end
