function bnd = prepare_mesh_segmentation(cfg, mri)

% PREPARE_MESH_SEGMENTATION
%
% See also PREPARE_MESH_MANUAL, PREPARE_MESH_HEADSHAPE, PREPARE_MESH_HEXAHEDRAL

% Copyrights (C) 2009, Robert Oostenveld
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
cfg.spmversion    = ft_getopt(cfg, 'spmversion', 'spm8');
cfg.method        = ft_getopt(cfg, 'method', 'projectmesh');
if all(isfield(mri, {'gray', 'white', 'csf'}))
  cfg.tissue      = ft_getopt(cfg, 'tissue', 'brain');    % set the default
  cfg.numvertices = ft_getopt(cfg, 'numvertices', 3000);  % set the default
else
  % do not set defaults for tissue and numvertices
  cfg.tissue      = ft_getopt(cfg, 'tissue');
  cfg.numvertices = ft_getopt(cfg, 'numvertices');
end

% check that SPM is on the path, try to add the preferred version
if strcmpi(cfg.spmversion, 'spm2'),
  ft_hastoolbox('SPM2',1);
elseif strcmpi(cfg.spmversion, 'spm8'),
  ft_hastoolbox('SPM8',1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try to determine the tissue (if not specified)

% special exceptional case first
if isempty(cfg.tissue) && numel(cfg.numvertices)==1 && isfield(mri,'white') && isfield(mri,'gray') && isfield(mri,'csf')
  mri=ft_datatype_segmentation(mri, 'segmentationstyle', 'probabilistic', 'hasbrain', 'yes');
  cfg.tissue='brain';
end

if isempty(cfg.tissue)
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  fn = fieldnames(mri);
  for i=1:numel(fn)
    if numel(mri.(fn{i}))==prod(mri.dim) && isfield(mri, [fn{i},'label'])
      segfield=fn{i};
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
elseif numel(cfg.tissue)~=numel(cfg.numvertices)
  error('you should specify the number of vertices for each tissue type');
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

for i =1:numel(cfg.tissue)
  if iscell(cfg.tissue)
    % the code below assumes that it is a probabilistic representation
    % for example {'brain', 'skull', scalp'}
    try
      seg = mri.(fixname(cfg.tissue{i}));
    catch
      error('Please specify cfg.tissue to correspond to tissue types in the segmented MRI')
    end
    tissue = cfg.tissue{i};
  else
    % this assumes that it is an indexed representation
    % for example [3 2 1]
    seg = (mri.seg==cfg.tissue(i));
    if isfield(mri, 'seglabel')
      try
        tissue = mri.seglabel{cfg.tissue(i)};
      catch
        error('Please specify cfg.tissue to correspond to (the name or number of) tissue types in the segmented MRI')
      end
    else
      tissue = sprintf('tissue %d', i);
    end
  end
  
  if strcmp(cfg.method, 'isosurface')
    fprintf('triangulating the outer boundary of compartment %d (%s) with the isosurface method\n', i, tissue);
  else
    fprintf('triangulating the outer boundary of compartment %d (%s) with %d vertices\n', i, tissue, cfg.numvertices(i));
  end
  
  % in principle it is possible to do volumesmooth and volumethreshold, but
  % the user is expected to prepare his segmentation outside this function
  % seg = volumesmooth(seg, nan, nan);
  
  % ensure that the segmentation is binary and that there is a single contiguous region
  seg = volumethreshold(seg, 0.5, tissue);
  
  % the function that generates the mesh will fail if there is a hole in the middle
  seg = volumefillholes(seg);
  
  switch cfg.method
    case 'isosurface'
      [tri, pnt] = isosurface(seg, 0.5);
      pnt = pnt(:,[2 1 3]); % Mathworks isosurface indexes differently
      
    case 'iso2mesh'
      ft_hastoolbox('iso2mesh', 1);
      
      opt = [];
      opt.radbound = 3; % set the target surface mesh element bounding sphere be <3 pixels in radius
      opt.maxnode = cfg.numvertices(i);
      opt.maxsurf = 1;
      
      [pnt, tri] = v2s(seg, 1, opt, 'cgalsurf');
      tri = tri(:,1:3);
      
    case 'projectmesh'
      [mrix, mriy, mriz] = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));
      ori(1) = mean(mrix(seg(:)));
      ori(2) = mean(mriy(seg(:)));
      ori(3) = mean(mriz(seg(:)));
      
      [pnt, tri] = triangulate_seg(seg, cfg.numvertices(i), ori);
      
    otherwise
      error('unsupported method "%s"', cfg.method);
  end % case
  
  numvoxels(i) = sum(find(seg(:))); % the number of voxels in this tissue
  
  bnd(i).pnt = ft_warp_apply(mri.transform, pnt);
  bnd(i).tri = tri;
  bnd(i).unit = mri.unit;
  
  
end % for each tissue

if strcmp(cfg.method, 'iso2surf')
  % order outside in (trying to be smart here)
  [dum, order] = sort(numvoxels,'descend');
  bnd = bnd(order);
  % clean up the triangulated meshes
  bnd = decouplesurf(bnd);
  % put them back in the original order
  [dum, order] = sort(order);
  bnd = bnd(order);
end

end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bnd = decouplesurf(bnd)
for ii = 1:length(bnd)-1
  % Despite what the instructions for surfboolean says, surfaces should
  % be ordered from inside-out!!
  [newnode, newelem] = surfboolean(bnd(ii+1).pnt,bnd(ii+1).tri,'decouple',bnd(ii).pnt,bnd(ii).tri);
  bnd(ii+1).tri = newelem(newelem(:,4)==2,1:3) - size(bnd(ii+1).pnt,1);
  bnd(ii+1).pnt = newnode(newnode(:,4)==2,1:3);
end % for
end %function
