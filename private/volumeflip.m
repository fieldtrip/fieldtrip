function [volume, flipvecout] = volumeflip(volume, flipvecin)

% VOLUMEFLIP
%
% See also VOLUMEPERMUTE, ALIGN_IJK2XYZ, ALIGN_XYZ2IJK

if nargin<2
  flipvecin = 'auto';
end

% do a low-level check on the input data
if ~isfield(volume, 'transform'), ft_error('the input volume needs a transformation matrix'); end
if ~isfield(volume, 'dim'),       ft_error('the input volume needs a dim field');             end

isrighthanded = det(volume.transform(1:3,1:3))>0;

% check the requested behavior
if ischar(flipvecin)
  switch flipvecin
    case {'left' 'lefthanded'}
      if isrighthanded
        ft_warning('left-handed axes system is requested, performing flip');
        flipvecin = [1 0 0];
      else
        flipvecin = [0 0 0];
      end
    case {'right' 'righthanded'}
      if ~isrighthanded
        ft_warning('right-handed axes system is requested, performing flip');
        flipvecin = [1 0 0];
      else
        flipvecin = [0 0 0];
      end
    case 'auto'
    otherwise
      ft_error('unsupported option specified');
  end
end

% get the largest values on the main diagonal of the transformation matrix
[volume, P] = volumepermute(volume);

% define some variable locally
T   = volume.transform;
dim = volume.dim;

% determine which fields can be flipped
fn = parameterselection('all', volume);

% pre-allocate
flipvecout = false(1,3);

% do the flipping on the numeric data and transformation matrix
for m = 1:3
  if ischar(flipvecin)
    flipvecout(m) = T(m,m)<0;
  else
    flipvecout(m) = flipvecin(m);
  end

  if flipvecout(m)
    % get the reflection matrix
    flipmat = eye(4); flipmat(m,m) = -1; flipmat(m,4) = dim(m)+1;
    for k = 1:numel(fn)
      tmp = volume.(fn{k});
      volume.(fn{k}) = flipdim(tmp, m);
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

