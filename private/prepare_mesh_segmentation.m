function mesh = prepare_mesh_segmentation(cfg, mri)

% PREPARE_MESH_SEGMENTATION
%
% The following configuration options can be specified for the iso2mesh method
%   cfg.maxsurf     = 1 = only use the largest disjointed surface
%                     0 = use all surfaces for that levelset
%   cfg.radbound    = a scalar indicating the radius of the target surface
%                     mesh element bounding sphere
%
% See also PREPARE_MESH_MANUAL, PREPARE_MESH_HEADSHAPE, PREPARE_MESH_HEXAHEDRAL,
% PREPARE_MESH_TETRAHEDRAL

% Copyrights (C) 2009-2021, Robert Oostenveld
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
cfg.spmversion    = ft_getopt(cfg, 'spmversion', 'spm12');
cfg.method        = ft_getopt(cfg, 'method', 'projectmesh');
cfg.maxsurf       = ft_getopt(cfg, 'maxsurf', 1);
cfg.radbound      = ft_getopt(cfg, 'radbound', 3);
if all(isfield(mri, {'gray', 'white', 'csf'}))
  % the default in this case is to combine gray, white and csf into brain, and segment that with 3000 vertices
  cfg.tissue      = ft_getopt(cfg, 'tissue', 'brain');    % set the default
  cfg.numvertices = ft_getopt(cfg, 'numvertices', 3000);  % set the default
else
  % do not set defaults for tissue and numvertices
  cfg.tissue      = ft_getopt(cfg, 'tissue');
  cfg.numvertices = ft_getopt(cfg, 'numvertices');
end

% check that the preferred SPM version is on the path
ft_hastoolbox(cfg.spmversion, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to determine the tissue (if not specified)

% special exceptional case first
if isempty(cfg.tissue) && numel(cfg.numvertices)==1 && isfield(mri, 'gray') && isfield(mri, 'white') && isfield(mri, 'csf')
  % combine the gray, white and csf into a single brain segmentation
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic', 'hasbrain', 'yes');
  cfg.tissue = 'brain';
end

if isempty(cfg.tissue)
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  fn = fieldnames(mri);
  for i=1:numel(fn)
    if numel(mri.(fn{i}))==prod(mri.dim) && isfield(mri, [fn{i}, 'label'])
      segfield = fn{i};
    end
  end
  if isfield(mri, [segfield 'label'])
    cfg.tissue = mri.([segfield 'label']);
    cfg.tissue = cfg.tissue(~cellfun(@isempty, cfg.tissue));
    fprintf('making mesh for %s tissue\n', cfg.tissue{:});
  else
    cfg.tissue=setdiff(unique(mri.(segfield)(:)),0);
    fprintf('making mesh for tissue with segmentation value %d\n', cfg.tissue);
  end
end

if ischar(cfg.tissue)
  % it should either be something like {'brain', 'skull', 'scalp'}, or something like [1 2 3]
  cfg.tissue = {cfg.tissue};
end

if numel(cfg.tissue)>1 && numel(cfg.numvertices)==1
  % use the same number of vertices for each tissue
  cfg.numvertices = repmat(cfg.numvertices, size(cfg.tissue));
end

if iscell(cfg.tissue)
  % the code below assumes that it is a probabilistic representation
  if any(strcmp(cfg.tissue, 'brain'))
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic', 'hasbrain', 'yes');
  else
    mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic');
  end
else
  % the code below assumes that it is an indexed representation
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the mesh extraction

for i=1:numel(cfg.tissue)
  
  if iscell(cfg.tissue)
    % the code below assumes that it is a probabilistic representation, for example {'brain', 'skull', scalp'}
    try
      seg = mri.(fixname(cfg.tissue{i}));
    catch
      ft_error('Please specify cfg.tissue to correspond to tissue types in the segmented MRI')
    end
    seglabel = cfg.tissue{i};
  elseif isnumeric(cfg.tissue) && isfield(mri, 'seg')
    % for backward compatibility with hard-coded seg field
    % this assumes that it is an indexed representation, for example [3 2 1]
    seg = (mri.seg==cfg.tissue(i));
    if isfield(mri, 'seglabel')
      try
        seglabel = mri.seglabel{cfg.tissue(i)};
      catch
        ft_error('specify cfg.tissue as the tissue name or tissue number')
      end
    else
      seglabel = sprintf('tissue %d', i);
    end
  elseif isnumeric(cfg.tissue) && isfield(mri, 'tissue')
    % for backward compatibility with hard-coded tissue field
    % this assumes that it is an indexed representation, for example [3 2 1]
    seg = (mri.tissue==cfg.tissue(i));
    if isfield(mri, 'tissuelabel')
      try
        seglabel = mri.tissuelabel{cfg.tissue(i)};
      catch
        ft_error('specify cfg.tissue as the tissue name or tissue number')
      end
    else
      seglabel = sprintf('tissue %d', i);
    end
  end
  
  if strcmp(cfg.method, 'isosurface')
    fprintf('triangulating the outer boundary of compartment %d (%s) with the isosurface method\n', i, seglabel);
  else
    fprintf('triangulating the outer boundary of compartment %d (%s) with %d vertices\n', i, seglabel, cfg.numvertices(i));
  end
  
  % in principle it is possible to do volumesmooth and volumethreshold, but
  % the user is expected to prepare his segmentation outside this function
  % seg = volumesmooth(seg, nan, nan);
  
  % ensure that the segmentation is binary and that there is a single contiguous region
  seg = volumethreshold(seg, 0.5, seglabel);
  
  % add a layer on all sides to ensure that the tissue can be meshed all the way up to the edges
  % this also ensures that the mesh at the bottom of the neck will be closed
  seg = volumepad(seg);
  
  % update the dimensions of the segmented volume to accomodate the padding
  dim = mri.dim+2;
  
  % update the homogenous transformation matrix to accomodate the shift that is due to the padding
  transform = mri.transform;
  shift = ft_warp_apply(transform, [1 1 1]) - ft_warp_apply(transform, [0 0 0]);
  transform(1,4) = transform(1,4) - shift(1);
  transform(2,4) = transform(2,4) - shift(2);
  transform(3,4) = transform(3,4) - shift(3);
  
  % the mesh generation will fail if there is a hole in the middle
  seg = volumefillholes(seg);
  
  switch cfg.method
    case 'isosurface'
      [tri, pos] = isosurface(seg);
      if ~isempty(cfg.numvertices)
        npos = cfg.numvertices(i);
        ntri = 2*(npos-2);
        [tri, pos] = reducepatch(tri, pos, ntri);
      end
      pos = pos(:,[2 1 3]); % Mathworks isosurface indexes differently
      
    case 'iso2mesh'
      % this requires the external iso2mesh toolbox
      ft_hastoolbox('iso2mesh', 1);
      
      opt = [];
      opt.radbound = cfg.radbound; % set the target surface mesh element bounding sphere be <3 pixels in radius
      opt.maxnode  = cfg.numvertices(i);
      opt.maxsurf  = cfg.maxsurf;
      
      method = 'cgalsurf';
      isovalues = 0.5;
      
      [pos, tri, regions, holes] = v2s(seg, isovalues, opt, method);
      pos = pos + 0.5; % undo the shift that has been introduced in vol2restrictedtri
      
      tri = tri(:,1:3);
      
    case 'projectmesh'
      [mrix, mriy, mriz] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
      ori(1) = mean(mrix(seg(:)));
      ori(2) = mean(mriy(seg(:)));
      ori(3) = mean(mriz(seg(:)));
      
      [pos, tri] = triangulate_seg(seg, cfg.numvertices(i), ori);
      
    otherwise
      ft_error('unsupported method "%s"', cfg.method);
  end % case
  
  numvoxels(i) = sum(find(seg(:))); % the number of voxels in this tissue
  
  mesh(i).pos  = ft_warp_apply(transform, pos);
  mesh(i).tri  = tri;
  mesh(i).unit = mri.unit;
  
end % for each tissue

if strcmp(cfg.method, 'iso2surf')
  % order outside in (trying to be smart here)
  [dum, order] = sort(numvoxels, 'descend');
  mesh = mesh(order);
  % clean up the triangulated meshes
  mesh = decouplesurf(mesh);
  % put them back in the original order
  [dum, order] = sort(order);
  mesh = mesh(order);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mesh = decouplesurf(mesh)
for ii = 1:length(mesh)-1
  % Despite what the instructions for surfboolean says, surfaces should be ordered from inside-out!!
  [newnode, newelem] = surfboolean(mesh(ii+1).pos, mesh(ii+1).tri, 'decouple', mesh(ii).pos,mesh(ii).tri);
  mesh(ii+1).tri = newelem(newelem(:,4)==2,1:3) - size(mesh(ii+1).pos,1);
  mesh(ii+1).pos = newnode(newnode(:,4)==2,1:3);
end % for
