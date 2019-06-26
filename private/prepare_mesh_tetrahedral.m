function mesh = prepare_mesh_tetrahedral(cfg, mri)

% PREPARE_MESH_TETRAHEDRAL
%
% See also PREPARE_MESH_MANUAL, PREPARE_MESH_HEADSHAPE,
% PREPARE_MESH_HEXAHEDRAL, PREPARE_MESH_SEGMENTATION

% Copyrights (C) 2016, Robert Oostenveld
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
% $Id$

% ensure that the input is consistent with what this function expects
mri = ft_checkdata(mri, 'datatype', {'volume', 'segmentation'}, 'hasunit', 'yes');

% get the default options
cfg.tissue = ft_getopt(cfg, 'tissue');

if isempty(cfg.tissue)
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  fn = fieldnames(mri);
  for i = 1:numel(fn)
    if (numel(mri.(fn{i})) == prod(mri.dim)) && (~strcmp(fn{i}, 'inside'))
      segfield = fn{i};
    end
  end
  cfg.tissue = setdiff(unique(mri.(segfield)(:)), 0);
end

if ischar(cfg.tissue)
  % it should either be something like {'brain', 'skull', 'scalp'}, or something like [1 2 3]
  cfg.tissue = {cfg.tissue};
end

if iscell(cfg.tissue)
  % the code below assumes that it is a probabilistic representation
  if any(strcmp(cfg.tissue, 'brain'))
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic', 'hasbrain', 'yes');
  else
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic');
  end
  % combine all tissue types
  seg = false(mri.dim);
  for i=1:numel(cfg.tissue)
    seg = seg | mri.(cfg.tissue{i});
  end
else
  % the code below assumes that it is an indexed representation
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  % combine all tissue types
  seg = (mri.seg>0);
end

% this requires the external iso2mesh toolbox
ft_hastoolbox('iso2mesh', 1);


[node, elem, face] = vol2mesh(seg, 1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3), 2, 2, 0,'cgalsurf');

mesh = keepfields(mri, {'coordsys', 'unit'});
mesh.pos = ft_warp_apply(mri.transform, node);
mesh.tet = elem(:,1:4);
