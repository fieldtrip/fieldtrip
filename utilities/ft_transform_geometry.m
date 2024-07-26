function [output] = ft_transform_geometry(transform, input, method)

% FT_TRANSFORM_GEOMETRY applies a homogeneous coordinate transformation to a
% structure with geometric information, for example a volume conduction model for the
% head, gradiometer of electrode structure containing EEG or MEG sensor positions and
% MEG coil orientations, a head shape or a source model.
%
% Use as
%   [output] = ft_transform_geometry(transform, input)
% where the transform should be a 4x4 homogeneous transformation matrix and the input
% data structure can be any of the FieldTrip data structures that describes
% geometrical data, or
%   [output] = ft_transform_geometry(transform, input, method)
% where the transform contains a set of parameters that can be converted into a 4x4 
% homogeneous transformation matrix, using one of the supported methods:
% 'rotate', 'scale', 'translate', 'rigidbody'. All methods require a 3-element vector
% as parameters, apart from rigidbody, which requires 6 parameters. 
%
% The units of the transformation matrix must be the same as the units in which the
% geometric object is expressed.
%
% The type of geometric object constrains the type of allowed transformations.
%
% For sensor arrays:
% If the input is an MEG gradiometer array, only a rigid-body translation plus
% rotation are allowed. If the input is an EEG electrode or fNIRS optodes array,
% global rescaling and individual axis rescaling is also allowed.
%
% For volume conduction models:
% If the input is a volume conductor model of the following type:
%   localspheres model
%   singleshell model with the spherical harmonic coefficients already computed
%   BEM model with system matrix already computed
%   FEM model with volumetric elements
% only a rigid-body translation plus rotation are allowed.
%
% If the input is a volume conductor model of the following type:
%   BEM model with the system matrix not yet computed
%   singleshell model with the spherical harmonic coefficients not yet computed
% rotation, translation, global rescaling and individual axis rescaling is allowed.
%
% If the input is a volume conductor model of the following type:
%   single sphere
%   concentric spheres
% rotation, translation and global rescaling is allowed.
%
% For source models, either defined as a 3D regular grid, a 2D mesh or unstructred
% point cloud, rotation, translation, global rescaling and individual axis rescaling
% is allowed.
%
% For anatomical MRIs and functional volumetric data, rotation, translation, global
% rescaling and individual axis rescaling are allowed.
%
% See also FT_WARP_APPLY, FT_HEADCOORDINATES, FT_SCALINGFACTOR

% Copyright (C) 2011-2024, Jan-Mathijs Schoffelen and Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_transform_geometry.m$

siz = size(transform);
if isequal(siz, [4 4])
  % this is OK
else
  % check whether the method has been specified, and to be consistent with
  % the input transform parameters, and create the transformation matrix
  if nargin<3
    ft_error('the first input argument is not a transformation matrix, hence a ''method'' should be specified');
  end
  switch method
    case {'scale' 'translate' 'rotate'}
      if numel(transform)~=3
        ft_error('the number of transformation parameters should be 3');
      end
    case 'rigidbody'
      if ~isequal(siz, [1 6]) && ~isequal(siz, [6 1])
        ft_error('the transformation parameters should contain six elements in a vector');
      end
    otherwise
      ft_error('unsupported method');
  end
  transform = feval(method, transform);
end

% determine the rotation matrix
rotation = eye(4);
rotation(1:3,1:3) = transform(1:3,1:3);

if any(abs(transform(4,:)-[0 0 0 1])>100*eps)
  ft_error('invalid transformation matrix');
end

% check whether the transformation includes a scaling operation
s = svd(rotation);
if abs(abs(det(rotation))-1)>1e-4
  % the transformation increases or decreases the overall volume, allow for
  % some numerical imprecision
  if any(abs(s./s(1)-1)>1e-3)
    % axes scale differently
    globalrescale = false;
    axesrescale   = true;
  else
    % global rescaling
    globalrescale = true;
    axesrescale   = true;
  end
else
  globalrescale = false;
  axesrescale   = false;
end

% check whether the input data combines well with the requested transformation
dtype = ft_datatype(input);
switch dtype
  case 'grad'
    if globalrescale || axesrescale, ft_error('only a rigid body transformation without rescaling is allowed'); end
  case {'mesh', 'mesh+label'}
    % there can be multiple meshes as a struct-array
    if isstruct(input) && length(input)>1
      for i=1:length(input)
        output(i) = ft_transform_geometry(transform, input(i));
      end
      return
    end
  otherwise
    % could be a volume conductor model with constrained transformation possibilities
    if ft_headmodeltype(input, 'singleshell') && isfield(input, 'forwpar') && (globalrescale || axesrescale)
      ft_error('only a rigid body transformation without rescaling is allowed');
    end
    if ft_headmodeltype(input, 'bem') && isfield(input, 'mat') && (globalrescale || axesrescale)
      ft_error('only a rigid body transformation without rescaling is allowed');
    end
    if ft_headmodeltype(input, 'localspheres') && (globalrescale || axesrescale)
      ft_error('only a rigid body transformation without rescaling is allowed');
    end
    if (isfield(input, 'tet') && isfield(input, 'stiff')) && (globalrescale || axesrescale)
      ft_error('only a rigid body transformation without rescaling is allowed');
    end
    if (isfield(input, 'hex') && isfield(input, 'stiff')) && (globalrescale || axesrescale)
      ft_error('only a rigid body transformation without rescaling is allowed');
    end
end

% tfields must be rotated, translated and scaled
% rfields must only be rotated
% mfields must be simply multiplied
% recfields must be recursed into
tfields   = {'pos' 'pnt' 'o' 'coilpos' 'elecpos' 'optopos' 'chanpos' 'chanposold' 'nas' 'lpa' 'rpa' 'zpoint'}; % apply rotation plus translation
rfields   = {'ori' 'nrm'     'coilori' 'elecori' 'optoori' 'chanori' 'chanoriold' 'mom' 'leadfield' 'filter'}; % only apply rotation
mfields   = {'transform'};                % plain matrix multiplication
recfields = {'fid' 'bnd' 'orig' 'dip'};   % recurse into these fields
% the field 'r' is not included here, because it applies to a volume
% conductor model, and scaling is not allowed, so r will not change.

fnames = fieldnames(input);
for k = 1:numel(fnames)
  % name = sprintf('%s.%s', inputname(2), fnames{k});
  if ~isempty(input.(fnames{k}))
    if any(strcmp(fnames{k}, tfields))
      % ft_info('applying transformation to %s', name);
      input.(fnames{k}) = apply(transform, input.(fnames{k}));
    elseif any(strcmp(fnames{k}, rfields))
      % ft_info('applying rotation to %s', name);
      input.(fnames{k}) = apply(rotation, input.(fnames{k}));
    elseif any(strcmp(fnames{k}, mfields))
      % ft_info('applying multiplication to %s', name);
      input.(fnames{k}) = transform*input.(fnames{k});
    elseif any(strcmp(fnames{k}, recfields))
      for j = 1:numel(input.(fnames{k}))
        % ft_info('recursing into %s', name);
        input.(fnames{k})(j) = ft_transform_geometry(transform, input.(fnames{k})(j));
      end
    else
      % do nothing
    end
  end
end
output = input;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that applies the homogeneous transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new] = apply(transform, old)

if iscell(old)
  % recurse into the cell-array
  new = cell(size(old));
  for i=1:numel(old)
    if ~isempty(old{i})
      new{i} = apply(transform, old{i});
    end
  end
  return;
end

[m, n] = size(old);
if m~=3 && n==3
  % each row is one position
  old(:,4) = 1;
  new = old * transform';
  new = new(:,1:3);
elseif m==3 && n~=3
  % each column is one position
  old(4,:) = 1;
  new = transform * old;
  new = new(1:3,:);
else
  % assume that each row is one position
  old(:,4) = 1;
  new = old * transform';
  new = new(:,1:3);
end
