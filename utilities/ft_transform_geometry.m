function [output] = ft_transform_geometry(transform, input)

% FT_TRANSFORM_GEOMETRY applies a homogeneous coordinate transformation to a
% structure with geometric information, for example a volume conduction model for the
% head, gradiometer of electrode structure containing EEG or MEG sensor positions and
% MEG coil orientations, a head shape or a source model.
%
% Use as
%   [output] = ft_transform_geometry(transform, input)
% where the transform should be a 4x4 homogenous transformation matrix and the
% input data structure can be any of the FieldTrip data structures that
% describes geometrical data.
%
% The units of the transformation matrix must be the same as the units in which the
% geometric object is expressed.
%
% Depending on the input object, the homogeneous transformation matrix should be
% limited to a rigid-body translation plus rotation (MEG-gradiometer array), or to a
% rigid-body translation plus rotation plus a global rescaling (volume conductor
% geometry).
%
% See also FT_WARP_APPLY, FT_HEADCOORDINATES, FT_SCALINGFACTOR

% Copyright (C) 2011, Jan-Mathijs Schoffelen
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

% determine the rotation matrix
rotation = eye(4);
rotation(1:3,1:3) = transform(1:3,1:3);

if any(abs(transform(4,:)-[0 0 0 1])>100*eps)
  ft_error('invalid transformation matrix');
end

% allow global rescaling for everything except MEG
if ft_senstype(input, 'meg')
  % allow for some numerical imprecision
  if (abs(det(rotation))-1)>1e-6
    ft_error('only a rigid body transformation without rescaling is allowed');
  end
end

% allow individual axis rescaling only for electrode positions and MRI-like geometries
if ~ft_senstype(input, 'eeg') && ~ft_datatype(input, 'volume')
  s = svd(transform(1:3,1:3));
  % allow for some numerical imprecision
  if any(abs(s./s(1)-1)>1e-3)
    ft_error('only a global rescaling is allowed');
  end
end

% tfields must be rotated, translated and scaled
% rfields must only be rotated
% mfields must be simply multipied
tfields   = {'pos' 'pnt' 'o' 'coilpos' 'elecpos' 'optopos' 'chanpos' 'chanposold' 'nas' 'lpa' 'rpa' 'zpoint'}; % apply rotation plus translation
rfields   = {'ori' 'nrm'     'coilori' 'elecori' 'optoori' 'chanori' 'chanoriold'                           }; % only apply rotation
mfields   = {'transform'};           % plain matrix multiplication
recfields = {'fid' 'bnd' 'orig'};    % recurse into these fields
% the field 'r' is not included here, because it applies to a volume
% conductor model, and scaling is not allowed, so r will not change.

fnames = fieldnames(input);
for k = 1:numel(fnames)
  if ~isempty(input.(fnames{k}))
    if any(strcmp(fnames{k}, tfields))
      input.(fnames{k}) = apply(transform, input.(fnames{k}));
    elseif any(strcmp(fnames{k}, rfields))
      input.(fnames{k}) = apply(rotation, input.(fnames{k}));
    elseif any(strcmp(fnames{k}, mfields))
      input.(fnames{k}) = transform*input.(fnames{k});
    elseif any(strcmp(fnames{k}, recfields))
      for j = 1:numel(input.(fnames{k}))
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
old(:,4) = 1;
new = old * transform';
new = new(:,1:3);
