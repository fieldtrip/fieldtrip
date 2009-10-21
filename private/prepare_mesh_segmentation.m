function bnd = prepare_mesh_segmentation(cfg, mri)

% PREPARE_MESH_SEGMENTATION
%
% See also PREPARE_MESH_MANUAL,PREPARE_MESH_HEADSHAPE

% Copyrights (C) 2009, Robert Oostenveld
%
% $Log: prepare_mesh_segmentation.m,v $
% Revision 1.4  2009/08/11 12:47:20  jansch
% fixed bug in assignment of default cfg.threshold
%
% Revision 1.3  2009/07/17 10:10:56  crimic
% added initial checks and variables consinstency
%
% Revision 1.2  2009/07/16 15:30:55  crimic
% fixed tiny error
%
% Revision 1.1  2009/07/13 14:45:06  crimic
% copy code of existing functions into stand-alone functions for inclusion in the prepare_mesh helper function
%

% some initial checks
cfg = checkconfig(cfg, 'forbidden', 'numcompartments');
if ~isfield(mri, 'tissue') && isfield(mri, 'gray'), cfg.tissue = 1; end
if ~isfield(cfg, 'threshold'), cfg.threshold = 0; end
  
fprintf('using the segmented MRI\n');

if ~isfield(mri, 'seg') && isequal(cfg.tissue, 1)
  mri.seg = zeros(size(mri.gray ));
  % construct the single image segmentation from the three probabilistic
  % tissue segmentations for csf, white and gray matter
  if isfield(mri, 'gray')
    fprintf('including gray matter in segmentation for brain compartment\n')
    mri.seg = mri.seg | (mri.gray>(cfg.threshold*max(mri.gray(:))));
  end
  if isfield(mri, 'white')
    fprintf('including white matter in segmentation for brain compartment\n')
    mri.seg = mri.seg | (mri.white>(cfg.threshold*max(mri.white(:))));
  end
  if isfield(mri, 'csf')
    fprintf('including CSF in segmentation for brain compartment\n')
    mri.seg = mri.seg | (mri.csf>(cfg.threshold*max(mri.csf(:))));
  end
  if ~strcmp(cfg.smooth, 'no'),
    % check whether the required SPM2 toolbox is available
    hastoolbox('spm2', 1);
    fprintf('smoothing the segmentation with a %d-pixel FWHM kernel\n',cfg.smooth);
    mri.seg = double(mri.seg);
    spm_smooth(mri.seg, mri.seg, cfg.smooth);
  end
  % threshold for the last time
  mri.seg = (mri.seg>(cfg.threshold*max(mri.seg(:))));
end

[mrix, mriy, mriz] = ndgrid(1:size(mri.seg,1), 1:size(mri.seg,2), 1:size(mri.seg,3));

% construct the triangulations of the boundaries from the segmented MRI
for i=1:length(cfg.tissue)
  fprintf('triangulating the boundary of compartment %d\n', i);
  seg = imfill((mri.seg==cfg.tissue(i)), 'holes');
  ori(1) = mean(mrix(seg(:)));
  ori(2) = mean(mriy(seg(:)));
  ori(3) = mean(mriz(seg(:)));
  [pnt, tri] = triangulate_seg(seg, cfg.numvertices(i), ori);
  % apply the coordinate transformation from voxel to head coordinates
  pnt(:,4) = 1;
  pnt = (mri.transform * (pnt'))';
  pnt = pnt(:,1:3);
  
  % convert the MRI surface points into the same units as the source/gradiometer
  scale = 1;
  switch cfg.sourceunits
    case 'mm'
      scale = scale * 1000;
    case 'cm'
      scale = scale * 100;
    case 'dm'
      scale = scale * 10;
    case 'm'
      scale = scale * 1;
    otherwise
      error('unknown physical dimension in cfg.sourceunits');
  end
  switch cfg.mriunits
    case 'mm'
      scale = scale / 1000;
    case 'cm'
      scale = scale / 100;
    case 'dm'
      scale = scale / 10;
    case 'm'
      scale = scale / 1;
    otherwise
      error('unknown physical dimension in cfg.mriunits');
  end
  if scale~=1
    fprintf('converting MRI surface points from %s into %s\n', cfg.sourceunits, cfg.mriunits);
    pnt = pnt* scale;
  end
  
  bnd(i).pnt = pnt;
  bnd(i).tri = tri;
  fprintf(['segmentation compartment %d of ' num2str(length(cfg.tissue)) ' completed\n'],i);
end

