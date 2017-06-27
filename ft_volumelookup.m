function [output] = ft_volumelookup(cfg, volume)

% FT_VOLUMELOOKUP can be used in to combine an anatomical or functional atlas with
% the source reconstruction results. You can use it for forward and reverse lookup.
%
% Given the ROI as anatomical or functional label, it looks up the locations and
% creates a mask (as a binary volume) based on the label. Given the ROI as point in
% the brain, it creates a sphere or box around that point. In these two case the
% function is to be used as:
%   mask = ft_volumelookup(cfg, volume)
%
% Given a binary volume that indicates a region of interest or a point of
% interest, it looks up the corresponding anatomical or functional labels
% from the atlas. In this case the function is to be used as:
%    labels = ft_volumelookup(cfg, volume)
%
% In both cases the input volume can be:
%   mri    is the output of FT_READ_MRI
%   source is the output of FT_SOURCEANALYSIS
%   stat   is the output of FT_SOURCESTATISTICS
%
% The configuration options for a mask according to an atlas:
%   cfg.inputcoord          = 'mni' or 'tal', coordinate system of the mri/source/stat
%   cfg.atlas               = string, filename of atlas to use, see FT_READ_ATLAS
%   cfg.roi                 = string or cell-array of strings, region(s) of interest from anatomical atlas
%
% The configuration options for a spherical/box mask around a point of interest:
%   cfg.roi                 = Nx3 vector, coordinates of the points of interest
%   cfg.sphere              = radius of each sphere in cm/mm dep on unit of input
%   cfg.box                 = Nx3 vector, size of each box in cm/mm dep on unit of input
%   cfg.round2nearestvoxel  = 'yes' or 'no' (default = 'no'), voxel closest to point of interest is calculated
%                             and box/sphere is centered around coordinates of that voxel
%
% The configuration options for labels from a mask:
%   cfg.inputcoord          = 'mni' or 'tal', coordinate system of the mri/source/stat
%   cfg.atlas               = string, filename of atlas to use, see FT_READ_ATLAS
%   cfg.maskparameter       = string, field in volume to be looked up, data in field should be logical
%   cfg.maxqueryrange       = number, should be 1, 3, 5 (default = 1)
%
% The configuration options for labels around a point of interest:
%   cfg.output              = 'label'
%   cfg.roi                 = Nx3 vector, coordinates of the points of interest
%   cfg.inputcoord          = 'mni' or 'tal', coordinate system of the mri/source/stat
%   cfg.atlas               = string, filename of atlas to use, see FT_READ_ATLAS
%   cfg.maxqueryrange       = number, should be 1, 3, 5 (default = 1)
%   cfg.round2nearestvoxel = 'yes' or 'no', voxel closest to point of interest is calculated (default = 'yes')
%
% The label output has a field "names", a field "count" and a field "usedqueryrange".
% To get a list of areas of the given mask you can do for instance:
%      [tmp ind] = sort(labels.count,1,'descend');
%      sel = find(tmp);
%      for j = 1:length(sel)
%        found_areas{j,1} = [num2str(labels.count(ind(j))) ': ' labels.name{ind(j)}];
%      end
% In the "found_areas" variable you can then see how many times which labels are
% found. Note that in the AFNI brick one location can have 2 labels.
%
% Dependent on the input coordinates and the coordinates of the atlas, the
% input MRI is transformed betweem MNI and Talairach-Tournoux coordinates
% See http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml for more details.
%
% See http://www.fieldtriptoolbox.org/template/atlas for a list of templates and
% atlasses that are included in the FieldTrip release.
%
% See also FT_READ_ATLAS, FT_SOURCEPLOT

% Copyright (C) 2008-2017, Robert Oostenveld, Ingrid Nieuwenhuis
% Copyright (C) 2013, Jan-Mathijs Schoffelen
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar volume
ft_preamble provenance volume
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the handling of the default cfg options is done further down
% the checking of the input data is done further down

cfg.maxqueryrange      = ft_getopt(cfg,'maxqueryrange', 1);
cfg.output             = ft_getopt(cfg,'output', []); % in future, cfg.output could be extended to support both 'label' and 'mask'

roi2mask   = 0;
mask2label = 0;
roi2label = 0;
if isfield(cfg, 'roi') && strcmp(cfg.output, 'label')
  roi2label = 1;
elseif isfield(cfg, 'roi')
  roi2mask = 1;
elseif isfield(cfg, 'maskparameter')
  mask2label = 1;
else
  error('you should either specify cfg.roi, or cfg.maskparameter')
end

if roi2mask
  % only for volume data
  volume = ft_checkdata(volume, 'datatype', 'volume');
  
  cfg.round2nearestvoxel = ft_getopt(cfg, 'round2nearestvoxel', 'no');
  
  isatlas = iscell(cfg.roi) || ischar(cfg.roi);
  ispoi   = isnumeric(cfg.roi);
  if isatlas+ispoi ~= 1
    error('do not understand cfg.roi')
  end
  
  if isatlas
    ft_checkconfig(cfg, 'forbidden', {'sphere', 'box'}, 'required',  {'atlas', 'inputcoord'});
  elseif ispoi
    ft_checkconfig(cfg, 'forbidden', {'atlas' 'inputcoord'});
    if isempty(ft_getopt(cfg, 'sphere')) && isempty(ft_getopt(cfg, 'box'))
      error('you should either specify cfg.sphere or cfg.box')
    end
  end
  
elseif mask2label || roi2label
  % convert to source representation (easier to work with)
  volume = ft_checkdata(volume, 'datatype', 'source');
  ft_checkconfig(cfg, 'required', {'atlas', 'inputcoord'});
  
  if isempty(intersect(cfg.maxqueryrange, [1 3 5]))
    error('incorrect query range, should be one of [1 3 5]');
  end
  
  if roi2label
    cfg.round2nearestvoxel = ft_getopt(cfg, 'round2nearestvoxel', 'yes');
  end
end


if roi2mask
  % determine location of each anatomical voxel in its own voxel coordinates
  dim = volume.dim;
  i = 1:dim(1);
  j = 1:dim(2);
  k = 1:dim(3);
  [I, J, K] = ndgrid(i, j, k);
  ijk = [I(:) J(:) K(:) ones(prod(dim),1)]';
  % determine location of each anatomical voxel in head coordinates
  xyz = volume.transform * ijk; % note that this is 4xN
  
  if isatlas
    if ischar(cfg.atlas)
      % assume it to represent a filename
      atlas = ft_read_atlas(cfg.atlas);
    else
      % assume cfg.atlas to be a struct, but it may have been converted
      % into a config object
      atlas = struct(cfg.atlas);
    end
    
    % determine which field(s) to use to look up the labels,
    % and whether these are boolean or indexed
    fn = fieldnames(atlas);
    isboolean = false(numel(fn),1);
    isindexed = false(numel(fn),1);
    for i=1:length(fn)
      if islogical(atlas.(fn{i})) && isequal(size(atlas.(fn{i})), atlas.dim)
        isboolean(i) = true;
      elseif isnumeric(atlas.(fn{i})) && isequal(size(atlas.(fn{i})), atlas.dim)
        isindexed(i) = true;
      end
    end
    if any(isindexed)
      % let the indexed prevail
      fn = fn(isindexed);
      isindexed = 1;
    elseif any(isboolean)
      % use the boolean
      fn = fn(isboolean);
      isindexed = 0;
    end
    
    if ischar(cfg.roi)
      cfg.roi = {cfg.roi};
    end
    
    if isindexed
      sel = zeros(0,2);
      for m = 1:length(fn)
        for i = 1:length(cfg.roi)
          tmp = find(strcmp(cfg.roi{i}, atlas.([fn{m},'label'])));
          sel = [sel; tmp m*ones(numel(tmp),1)];
        end
      end
      fprintf('found %d matching anatomical labels\n', size(sel,1));
      
      % this is to accommodate for multiple parcellations:
      % the brick refers to the parcellationname
      % the value refers to the value within the given parcellation
      brick = sel(:,2);
      value = sel(:,1);
      
      % convert between MNI head coordinates and TAL head coordinates
      % coordinates should be expressed compatible with the atlas
      if     strcmp(cfg.inputcoord, 'mni') && strcmp(atlas.coordsys, 'tal')
        xyz(1:3,:) = mni2tal(xyz(1:3,:));
      elseif strcmp(cfg.inputcoord, 'mni') && strcmp(atlas.coordsys, 'mni')
        % nothing to do
      elseif strcmp(cfg.inputcoord, 'tal') && strcmp(atlas.coordsys, 'tal')
        % nothing to do
      elseif strcmp(cfg.inputcoord, 'tal') && strcmp(atlas.coordsys, 'mni')
        xyz(1:3,:) = tal2mni(xyz(1:3,:));
      elseif ~strcmp(cfg.inputcoord, atlas.coordsys)
        error('there is a mismatch between the coordinate system in the atlas and the coordinate system in the data, which cannot be resolved');
      end
      
      % determine location of each anatomical voxel in atlas voxel coordinates
      ijk = atlas.transform \ xyz;
      ijk = round(ijk(1:3,:))';
      
      inside_vol = ijk(:,1)>=1 & ijk(:,1)<=atlas.dim(1) & ...
        ijk(:,2)>=1 & ijk(:,2)<=atlas.dim(2) & ...
        ijk(:,3)>=1 & ijk(:,3)<=atlas.dim(3);
      inside_vol = find(inside_vol);
      
      % convert the selection inside the atlas volume into linear indices
      ind = sub2ind(atlas.dim, ijk(inside_vol,1), ijk(inside_vol,2), ijk(inside_vol,3));
      
      brick_val = cell(1,numel(brick));
      % search the bricks for the value of each voxel
      for i=1:numel(brick_val)
        brick_val{i} = zeros(prod(dim),1);
        brick_val{i}(inside_vol) = atlas.(fn{brick(i)})(ind);
      end
      
      mask = zeros(prod(dim),1);
      for i=1:numel(brick_val)
        %fprintf('constructing mask for %s\n', atlas.descr.name{sel(i)});
        mask = mask | (brick_val{i}==value(i));
      end
    else
      error('support for atlases that have a probabilistic segmentationstyle is not supported yet');
      % NOTE: this may be very straightforward indeed: the mask is just the
      % logical or of the specified rois.
    end
    
  elseif ispoi
    
    if istrue(cfg.round2nearestvoxel)
      for i=1:size(cfg.roi,1)
        cfg.roi(i,:) = poi2voi(cfg.roi(i,:), xyz);
      end
    end
    
    % sphere(s)
    if isfield(cfg, 'sphere')
      mask = zeros(1,prod(dim));
      for i=1:size(cfg.roi,1)
        dist = sqrt( (xyz(1,:) - cfg.roi(i,1)).^2 + (xyz(2,:) - cfg.roi(i,2)).^2 + (xyz(3,:) - cfg.roi(i,3)).^2 );
        mask = mask | (dist <= cfg.sphere(i));
      end
      % box(es)
    elseif isfield(cfg, 'box')
      mask = zeros(1, prod(dim));
      for i=1:size(cfg.roi,1)
        mask = mask | ...
          (xyz(1,:) <= (cfg.roi(i,1) + cfg.box(i,1)./2) & xyz(1,:) >= (cfg.roi(i,1) - cfg.box(i,1)./2)) & ...
          (xyz(2,:) <= (cfg.roi(i,2) + cfg.box(i,2)./2) & xyz(2,:) >= (cfg.roi(i,2) - cfg.box(i,2)./2)) & ...
          (xyz(3,:) <= (cfg.roi(i,3) + cfg.box(i,3)./2) & xyz(3,:) >= (cfg.roi(i,3) - cfg.box(i,3)./2));
      end
    end
  end
  
  mask = reshape(mask, dim);
  fprintf('%i voxels in mask, which is %.3f %% of total volume\n', sum(mask(:)), 100*mean(mask(:)));
  output = mask;
  
elseif mask2label || roi2label
  if ischar(cfg.atlas)
    % assume it to represent a filename
    atlas = ft_read_atlas(cfg.atlas);
  else
    % assume cfg.atlas to be a struct, but it may have been converted
    % into a config object
    atlas = struct(cfg.atlas);
  end
  
  % determine which field(s) to use to look up the labels,
  % and whether these are boolean or indexed
  fn = fieldnames(atlas);
  isboolean = false(numel(fn),1);
  isindexed = false(numel(fn),1);
  for i=1:length(fn)
    if islogical(atlas.(fn{i})) && isequal(size(atlas.(fn{i})), atlas.dim)
      isboolean(i) = true;
    elseif isnumeric(atlas.(fn{i})) && isequal(size(atlas.(fn{i})), atlas.dim)
      isindexed(i) = true;
    end
  end
  if any(isindexed)
    % let the indexed prevail
    fn = fn(isindexed);
    isindexed = 1;
  elseif any(isboolean)
    % use the boolean
    fn = fn(isboolean);
    isindexed = 0;
  end
  
  labels.name = cell(0,1);
  for k = 1:numel(fn)
    % ensure that they are concatenated as column
    tmp = atlas.([fn{k},'label']);
    labels.name = cat(1, labels.name(:), tmp(:));
  end
  labels.name{end+1} = 'no_label_found';
  labels.count = zeros(length(labels.name),1);
  for iLab = 1:length(labels.name)
    labels.usedqueryrange{iLab} = [];
  end
  
  if mask2label
    sel = find(volume.(cfg.maskparameter)(:));
  elseif roi2label
    if istrue(cfg.round2nearestvoxel)
      % determine location of each anatomical voxel in head coordinates
      xyz = [volume.pos ones(size(volume.pos,1),1)]'; % note that this is 4xN
      for i=1:size(cfg.roi,1)
        cfg.roi(i,:) = poi2voi(cfg.roi(i,:), xyz);
      end
    end % round2nearestvoxel
    sel = find(ismember(volume.pos, cfg.roi, 'rows')==1);
  end
  for iVox = 1:length(sel)
    usedQR = 1;
    label = atlas_lookup(atlas, [volume.pos(sel(iVox),1) volume.pos(sel(iVox),2) volume.pos(sel(iVox),3)], 'inputcoord', cfg.inputcoord, 'queryrange', 1);
    if isempty(label) && cfg.maxqueryrange > 1
      label = atlas_lookup(atlas, [volume.pos(sel(iVox),1) volume.pos(sel(iVox),2) volume.pos(sel(iVox),3)], 'inputcoord', cfg.inputcoord, 'queryrange', 3);
      usedQR = 3;
    end
    if isempty(label) && cfg.maxqueryrange > 3
      label = atlas_lookup(atlas, [volume.pos(sel(iVox),1) volume.pos(sel(iVox),2) volume.pos(sel(iVox),3)], 'inputcoord', cfg.inputcoord, 'queryrange', 5);
      usedQR = 5;
    end
    if isempty(label)
      label = {'no_label_found'};
    elseif length(label) == 1
      label = {label};
    end
    
    ind_lab = [];
    for iLab = 1:length(label)
      ind_lab = [ind_lab find(strcmp(label{iLab}, labels.name))];
    end
    
    labels.count(ind_lab) = labels.count(ind_lab) + (1/length(ind_lab));
    for iFoundLab = 1:length(ind_lab)
      if isempty(labels.usedqueryrange{ind_lab(iFoundLab)})
        labels.usedqueryrange{ind_lab(iFoundLab)} = usedQR;
      else
        labels.usedqueryrange{ind_lab(iFoundLab)} = [labels.usedqueryrange{ind_lab(iFoundLab)} usedQR];
      end
    end
  end %iVox
  
  output = labels;
  
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION point of interest to voxel of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function voi = poi2voi(poi, xyz)
xmin = min(abs(xyz(1,:) - poi(1))); xcl = round(abs(xyz(1,:) - poi(1))) == round(xmin);
ymin = min(abs(xyz(2,:) - poi(2))); ycl = round(abs(xyz(2,:) - poi(2))) == round(ymin);
zmin = min(abs(xyz(3,:) - poi(3))); zcl = round(abs(xyz(3,:) - poi(3))) == round(zmin);
xyzcls = xcl + ycl + zcl; ind_voi = xyzcls == 3;
if sum(ind_voi) > 1
  fprintf('%i voxels at same distance of poi, taking first voxel\n', sum(ind_voi))
  ind_voi_temp = find(ind_voi); ind_voi_temp = ind_voi_temp(1);
  ind_voi = zeros(size(ind_voi));
  ind_voi(ind_voi_temp) = 1;
  ind_voi = logical(ind_voi);
end
voi = xyz(1:3,ind_voi);
fprintf('coordinates of voi: %.1f  %.1f %.1f\n', voi(1), voi(2), voi(3));
