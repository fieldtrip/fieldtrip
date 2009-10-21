function [output] = volumelookup(cfg, volume)

% VOLUMELOOKUP can be used in two directions. (1) It creates a mask
% according to a label from a given atlas, or a sphere or box around a
% point of interest. (2) It gives the labels from a given atlas according
% to a given mask.
%
% Use as
%   (1) [mask]  = volumelookup(cfg, volume)
% or
%   (2) [labels] = volumelookup(cfg, volume)
%
% where volume can be:
%   mri    is the output of READ_FCDC_MRI
%   source is the output of SOURCEANALYSIS
%   stat   is the output of SOURCESTATISTICS
%
% configuration options for a mask according to an atlas:
%   cfg.inputcoord = 'mni' or 'tal', coordinate system of the mri/source/stat
%   cfg.atlas      = string, filename of atlas to use, either the AFNI
%                     brik file that is available from http://afni.nimh.nih.gov/afni/doc/misc/ttatlas_tlrc,
%                     or the WFU atlasses available from  http://fmri.wfubmc.edu. see PREPARE_ATLAS
%   cfg.roi        = string or cell of strings, region(s) of interest from anatomical atlas
% configuration options for a spherical/box mask around a point of interest:
%   cfg.roi                = Nx3 vector, coordinates of the points of interest
%   cfg.sphere             = radius of each sphere in cm/mm dep on unit of input
%   cfg.box                = Nx3 vector, size of each box in cm/mm dep on unit of input
%   cfg.round2nearestvoxel = 'yes' or 'no' (default = 'no'), voxel closest to point of interest is calculated
%                             and box/sphere is centered around coordinates of that voxel
% configuration options for labels from a mask:
%   cfg.inputcoord    = 'mni' or 'tal', coordinate system of the mri/source/stat
%   cfg.atlas         = string, filename of atlas to use, either the AFNI
%                        brik file that is available from http://afni.nimh.nih.gov/afni/doc/misc/ttatlas_tlrc,
%                        or the WFU atlasses available from http://fmri.wfubmc.edu. see PREPARE_ATLAS
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

% Copyright (C) 2008, Robert Oostenveld, Ingrid Nieuwenhuis
%
% $Log: volumelookup.m,v $
% Revision 1.5  2009/08/03 15:20:19  ingnie
% updated help
%
% Revision 1.4  2009/01/29 13:54:57  ingnie
% changed help
%
% Revision 1.3  2009/01/27 10:09:49  ingnie
% added mask2label functionality
%
% Revision 1.2  2008/12/16 20:42:33  roboos
% converted from dos to unix text file
%
% Revision 1.1  2008/12/05 13:39:38  ingnie
% new implementation based om atlas_mask, added possibility to make spherical and box masks.
%

fieldtripdefs

roi2mask = 0;
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
  volume = checkdata(volume, 'datatype', 'volume');
  
  % set the defaults
  if ~isfield(cfg, 'round2nearestvoxel'),  cfg.round2nearestvoxel = 'no';  end

  if iscell(cfg.roi) || ischar(cfg.roi)
    checkconfig(cfg, 'forbidden', {'sphere' 'box'}, 'required', {'atlas' 'inputcoord'});
    isatlas = 1;
    ispoi = 0;
  elseif isnumeric(cfg.roi)
    checkconfig(cfg, 'forbidden', {'atlas' 'inputcoord'});
    isatlas = 0;
    ispoi = 1;
  else
    error('do not understand cfg.roi')
  end

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
    atlas = prepare_atlas(cfg.atlas);

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

    if strcmp(cfg.round2nearestvoxel, 'yes')
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
    else
      error('either specify cfg.sphere or cfg.box')
    end
  end

  mask = reshape(mask, dim);
  fprintf('%i voxels in mask, which is %.3f %% of total volume\n', sum(mask(:)), 100*mean(mask(:)));
  output = mask;
  
elseif mask2label
  % convert to source representation (easier to work with)
  volume = checkdata(volume, 'datatype', 'source');
  checkconfig(cfg, 'required', {'atlas' 'inputcoord'});

  % set defaults
  if ~isfield(cfg, 'maxqueryrange'),  cfg.maxqueryrange = 1;  end

  if isempty(intersect(cfg.maxqueryrange, [1 3 5]))
    error('incorrect query range, should be one of [1 3 5]');
  end

  atlas = prepare_atlas(cfg.atlas);
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

