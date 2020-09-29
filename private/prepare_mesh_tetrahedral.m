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
  % this requires an indexed segmentation, in which case cfg.tissue will be
  % a numeric vector, that will be defined below
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
end

fn = fieldnames(mri);
fn(strcmp(fn, 'inside')) = [];
[indexed, probabilistic] = determine_segmentationstyle(mri, fn, mri.dim);
if sum(indexed)==0
  % what to do?
elseif sum(indexed)==1
  segfield = fn{indexed};
elseif sum(indexed)>1
  ft_error('The input segmentation contains more than indexed segmented volume. Unsupported');
end

if isempty(cfg.tissue)
  cfg.tissue = setdiff(unique(mri.(segfield)(:)), 0);
end

if ischar(cfg.tissue)
  % it should either be something like {'brain', 'skull', 'scalp'}, or something like [1 2 3]
  cfg.tissue = {cfg.tissue};
end

% from here, cfg.tissue is either a cell-array of strings, or a numeric vector
if iscell(cfg.tissue)
  % the code below assumes that it is a probabilistic representation
  if any(strcmp(cfg.tissue, 'brain'))
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic', 'hasbrain', 'yes');
  else
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic');
  end
  % combine all tissue types
  seg      = false(mri.dim);
  seglabel = cell(numel(cfg.tissue),1); 
  for i=1:numel(cfg.tissue)
    seg = seg + double(mri.(cfg.tissue{i})).*i;
    seglabel{i} = cfg.tissue{i};
  end
else
  % the code below assumes that it is an indexed representation and that
  % cfg.tissue is numeric. FIXME: it also assumes that the indexed volume
  % is called 'seg', and that the label is called 'seglabel'. This is
  % probably a remnant of old code.
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  % keep the tissue types as they are
  seg      = mri.seg;
  seglabel = mri.seglabel;
end

% this requires the external iso2mesh toolbox
ft_hastoolbox('iso2mesh', 1);

[node, elem, face] = vol2mesh(uint8(seg), 1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3), 1, 1, 0, 'cgalmesh');
elem(:,1:4)        = meshreorient(node(:,1:3), elem(:,1:4));

mesh = keepfields(mri, {'coordsys', 'unit'});
mesh.pos = ft_warp_apply(mri.transform, node(:,1:3)+1); % offset of 1 is needed because indexing is 0-based?
mesh.tet = elem(:,1:4);
mesh.tissue = elem(:,5);
if exist('seglabel', 'var')
  mesh.tissuelabel = seglabel;
end