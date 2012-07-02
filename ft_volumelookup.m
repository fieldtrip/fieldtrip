function [output] = ft_volumelookup(cfg, volume)

% FT_VOLUMELOOKUP can be used in to combine an anatomical or functional
% atlas with source recunstruction. You can use it for forward and reverse
% lookup.
%
% Given the anatomical or functional label, it looks up the locations and
% creates a mask (as a binary volume) based on the label, or creates a
% sphere or box around a point of interest. In this case the function is to
% be used as:
%   mask = ft_volumelookup(cfg, volume)
%
% Given a binary volume that indicates a region of interest, it looks up
% the corresponding anatomical or functional labels from a given atlas. In
% this case the function is to be used as follows:
%    labels = ft_volumelookup(cfg, volume)
%
% In both cases the input volume can be:
%   mri    is the output of FT_READ_MRI
%   source is the output of FT_SOURCEANALYSIS
%   stat   is the output of FT_SOURCESTATISTICS
%
% The configuration options for a mask according to an atlas:
%   cfg.inputcoord = 'mni' or 'tal', coordinate system of the mri/source/stat
%   cfg.atlas      = string, filename of atlas to use, either the AFNI
%                     brik file that is available from http://afni.nimh.nih.gov/afni/doc/misc/ttatlas_tlrc,
%                     or the WFU atlasses available from  http://fmri.wfubmc.edu. see FT_PREPARE_ATLAS
%   cfg.roi        = string or cell of strings, region(s) of interest from anatomical atlas
%
% The configuration options for a spherical/box mask around a point of interest:
%   cfg.roi                = Nx3 vector, coordinates of the points of interest
%   cfg.sphere             = radius of each sphere in cm/mm dep on unit of input
%   cfg.box                = Nx3 vector, size of each box in cm/mm dep on unit of input
%   cfg.round2nearestvoxel = 'yes' or 'no' (default = 'no'), voxel closest to point of interest is calculated
%                             and box/sphere is centered around coordinates of that voxel
%
% The configuration options for labels from a mask:
%   cfg.inputcoord    = 'mni' or 'tal', coordinate system of the mri/source/stat
%   cfg.atlas         = string, filename of atlas to use, either the AFNI
%                        brik file that is available from http://afni.nimh.nih.gov/afni/doc/misc/ttatlas_tlrc,
%                        or the WFU atlasses available from http://fmri.wfubmc.edu. see FT_PREPARE_ATLAS
%   cfg.maskparameter = string, field in volume to be lookedup, data in field should be logical
%   cfg.maxqueryrange = number, should be 1, 3, 5 (default = 1)
%
% The label output has a field "names", a field "count" and a field "usedqueryrange"
% To get a list of areas of the given mask you can do for instance:
%      [tmp ind] = sort(labels.count,1,'descend');
%      sel = find(tmp);
%      for j = 1:length(sel)
%        found_areas{j,1} = [num2str(labels.count(ind(j))) ': ' labels.name{ind(j)}];
%      end
% in found_areas you can then see how many times which labels are found
% NB in the AFNI brick one location can have 2 labels!
%
% Dependent on the input coordinates and the coordinates of the atlas, the
% input MRI is transformed betweem MNI and Talairach-Tournoux coordinates
% See http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml for more details.
%
% See also FT_PREPARE_ATLAS, FT_SOURCEPLOT

% Copyright (C) 2008, Robert Oostenveld, Ingrid Nieuwenhuis
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar volume

% the handling of the default cfg options is done further down
% the checking of the input data is done further down

cfg.maxqueryrange      = ft_getopt(cfg,'maxqueryrange', 1);

roi2mask   = 0;
mask2label = 0;
if isfield(cfg, 'roi');
  roi2mask = 1;
elseif isfield(cfg, 'maskparameter')
  mask2label = 1;
else
  error('either specify cfg.roi, or cfg.maskparameter')
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
    ft_checkconfig(cfg, 'forbidden', {'sphere' 'box'}, ...
                        'required',  {'atlas' 'inputcoord'});  
  elseif ispoi
    ft_checkconfig(cfg, 'forbidden', {'atlas' 'inputcoord'});
    if isempty(ft_getopt(cfg, 'sphere')) && isempty(ft_getopt(cfg, 'box'))
      % either needs to be there
      error('either specify cfg.sphere or cfg.box')
    end  
  end
  
elseif mask2label
  % convert to source representation (easier to work with)
  volume = ft_checkdata(volume, 'datatype', 'source');
  ft_checkconfig(cfg, 'required', {'atlas' 'inputcoord'});
  
  if isempty(intersect(cfg.maxqueryrange, [1 3 5]))
    error('incorrect query range, should be one of [1 3 5]');
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
    atlas = ft_prepare_atlas(cfg);

    if ischar(cfg.roi)
      cfg.roi = {cfg.roi};
    end

    sel = [];
    for i = 1:length(cfg.roi)
      sel = [sel; strmatch(cfg.roi{i}, atlas.descr.name, 'exact')];
    end

    fprintf('found %d matching anatomical labels\n', length(sel));

    brick = atlas.descr.brick(sel);
    value = atlas.descr.value(sel);

    % convert between MNI head coordinates and TAL head coordinates
    % coordinates should be expressed compatible with the atlas
    if     strcmp(cfg.inputcoord, 'mni') && strcmp(atlas.coord, 'tal')
      xyz(1:3,:) = mni2tal(xyz(1:3,:));
    elseif strcmp(cfg.inputcoord, 'mni') && strcmp(atlas.coord, 'mni')
      % nothing to do
    elseif strcmp(cfg.inputcoord, 'tal') && strcmp(atlas.coord, 'tal')
      % nothing to do
    elseif strcmp(cfg.inputcoord, 'tal') && strcmp(atlas.coord, 'mni')
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

elseif mask2label
  atlas = ft_prepare_atlas(cfg.atlas);
  sel = find(volume.(cfg.maskparameter)(:));
  labels.name = atlas.descr.name;
  labels.name{end+1} = 'no_label_found';
  labels.count = zeros(length(labels.name),1);
  for iLab = 1:length(labels.name)
    labels.usedqueryrange{iLab} = [];
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
ft_postamble trackconfig
ft_postamble callinfo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION point of interest to voxel of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function voi = poi2voi(poi, xyz)
xmin = min(abs(xyz(1,:) - poi(1))); xcl = round(abs(xyz(1,:) - poi(1))) == round(xmin);
ymin = min(abs(xyz(2,:) - poi(2))); ycl = round(abs(xyz(2,:) - poi(2))) == round(ymin);
zmin = min(abs(xyz(3,:) - poi(3))); zcl = round(abs(xyz(3,:) - poi(3))) == round(zmin);
xyzcls = xcl + ycl + zcl; ind_voi = xyzcls == 3;
if sum(ind_voi) > 1;
  fprintf('%i voxels at same distance of poi, taking first voxel\n', sum(ind_voi))
  ind_voi_temp = find(ind_voi); ind_voi_temp = ind_voi_temp(1);
  ind_voi = zeros(size(ind_voi));
  ind_voi(ind_voi_temp) = 1;
  ind_voi = logical(ind_voi);
end
voi = xyz(1:3,ind_voi);
fprintf('coordinates of voi: %.1f  %.1f %.1f\n', voi(1), voi(2), voi(3));

