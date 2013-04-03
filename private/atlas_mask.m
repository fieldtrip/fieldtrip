function [mask] = atlas_mask(atlas, mri, label, varargin)

% ATLAS_MASK creates a mask that can be used in visualizing a functional
% and/or anatomical MRI volume.
%
% Use as
%   atlas = atlas_init;
%   mask  = atlas_mask(atlas, mri, label, ...);
%
% Optinal input arguments should come in key-value pairs and can include
%   'inputcoord'   = 'mni' or 'tal' (default = []);
%
% Dependent on the input coordinates and the coordinates of the atlas, the
% input MRI is transformed betweem MNI and Talairach-Tournoux coordinates
% See http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml for more details.
%
% See also ATLAS_INIT, ATLAS_LOOKUP

% Copyright (C) 2005-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
% $Id$

% get the optional input arguments
inputcoord = ft_getopt(varargin, 'inputcoord');

if isempty(inputcoord)
  error('you must specify inputcoord');
end

if ischar(label)
  label = {label};
end

sel = [];
for i=1:length(label)
  sel = [sel; find(strcmp(label{i}, atlas.descr.name))];
end

fprintf('found %d matching anatomical labels\n', length(sel));

brick = atlas.descr.brick(sel);
value = atlas.descr.value(sel);

if isfield(mri, 'transform') && isfield(mri, 'dim')
  dim = mri.dim;
  % determine location of each anatomical voxel in its own voxel coordinates
  i = 1:dim(1);
  j = 1:dim(2);
  k = 1:dim(3);
  [I, J, K] = ndgrid(i, j, k);
  ijk = [I(:) J(:) K(:) ones(prod(dim),1)]';
  % determine location of each anatomical voxel in head coordinates
  xyz = mri.transform * ijk; % note that this is 4xN
elseif isfield(mri, 'pos')
  % the individual positions of every grid point are specified
  npos = size(mri.pos,1);
  dim  = [npos 1];
  xyz  = [mri.pos ones(npos,1)]';  % note that this is 4xN
else
  error('could not determine whether the input describes a volume or a source');
end

% convert between MNI head coordinates and TAL head coordinates
% coordinates should be expressed compatible with the atlas
if     strcmp(inputcoord, 'mni') && strcmp(atlas.coord, 'tal')
  xyz(1:3,:) = mni2tal(xyz(1:3,:));
elseif strcmp(inputcoord, 'mni') && strcmp(atlas.coord, 'mni')
  % nothing to do
elseif strcmp(inputcoord, 'tal') && strcmp(atlas.coord, 'tal')
  % nothing to do
elseif strcmp(inputcoord, 'tal') && strcmp(atlas.coord, 'mni')
  xyz(1:3,:) = tal2mni(xyz(1:3,:));
end

% determine location of each anatomical voxel in atlas voxel coordinates
ijk = inv(atlas.transform) * xyz;
ijk = round(ijk(1:3,:))';

inside_vol = ijk(:,1)>=1 & ijk(:,1)<=atlas.dim(1) & ...
  ijk(:,2)>=1 & ijk(:,2)<=atlas.dim(2) & ...
  ijk(:,3)>=1 & ijk(:,3)<=atlas.dim(3);
inside_vol = find(inside_vol);

% convert the selection inside the atlas volume into linear indices
ind = sub2ind(atlas.dim, ijk(inside_vol,1), ijk(inside_vol,2), ijk(inside_vol,3));

brick0_val = zeros(prod(dim),1);
brick1_val = zeros(prod(dim),1);
% search the two bricks for the value of each voxel
brick0_val(inside_vol) = atlas.brick0(ind);
brick1_val(inside_vol) = atlas.brick1(ind);

mask = zeros(prod(dim),1);
for i=1:length(sel)
  fprintf('constructing mask for %s\n', atlas.descr.name{sel(i)});
  if brick(i)==0
    mask = mask | (brick0_val==value(i));
  elseif brick(i)==1
    mask = mask | (brick1_val==value(i));
  end
end
mask = reshape(mask, dim);

fprintf('masked %.1f %% of total volume\n', 100*mean(mask(:)));

