function [volume, flipvec] = volumeflip(volume, flipvecin)

if nargin<2
  flipvecin = 'auto';
end

% do a low-level check on the input data
if ~isfield(volume, 'transform'), error('the input volume needs a transformation matrix'); end
if ~isfield(volume, 'dim'),       error('the input volume needs a dim field');             end

if isfield(volume, 'outside')
end

isrighthanded = det(volume.transform(1:3,1:3))>0;

% check the requested behavior
if ischar(flipvecin)
  switch flipvecin
    case {'left' 'lefthanded'}
      if isrighthanded
        warning('left-handed axes system is requested, performing flip');
        flipvecin = [1 0 0];
      else
        flipvecin = [0 0 0];
      end
    case {'right' 'righthanded'}
      if ~isrighthanded
        warning('right-handed axes system is requested, performing flip');
        flipvecin = [1 0 0];
      else
        flipvecin = [0 0 0];
      end
    case 'auto'
    otherwise
      error('unsupported option specified');
  end
end

% get the largest values on the main diagonal of the transformation matrix
[volume, P] = volumepermute(volume);

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

% pre-allocate
flipvec = false(1,3);

% do the flipping on the numeric data and transformation matrix
for m = 1:3
  if ischar(flipvecin)
    flipvec(m) = T(m,m)<0;
  else
    flipvec(m) = flipvecin(m);
  end
  
  if flipvec(m)
    % get the reflection matrix
    flipmat = eye(4); flipmat(m,m) = -1; flipmat(m,4) = dim(m)+1; 
    for k = 1:numel(fnames)
      volume = setsubfield(volume, fnames{k}, flipdim(getsubfield(volume, fnames{k}), m));
    end
    T = T*flipmat;
  else
    % do nothing
  end
end

% update the transformation matrix
volume.transform = T;

% permute back
volume = volumepermute(volume, P);

