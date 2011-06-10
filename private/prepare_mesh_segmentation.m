function bnd = prepare_mesh_segmentation(cfg, mri)

% PREPARE_MESH_SEGMENTATION
%
% See also PREPARE_MESH_MANUAL,PREPARE_MESH_HEADSHAPE

% Copyrights (C) 2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

if ~isfield(cfg, 'spmversion'), cfg.spmversion = 'spm8'; end

% smooth functional parameters, excluding anatomy and inside
if isfield(cfg, 'smooth') && ~strcmp(cfg.smooth, 'no'),
  % check that SPM is on the path, try to add the preferred version
  if strcmpi(cfg.spmversion, 'spm2'),
    ft_hastoolbox('SPM2',1);
  elseif strcmpi(cfg.spmversion, 'spm8'),
    ft_hastoolbox('SPM8',1);
  end
end

% some initial checks
cfg = ft_checkconfig(cfg, 'forbidden', 'numcompartments');
if ~isfield(mri, 'tissue') && any(ismember(fieldnames(mri), {'gray' 'brain' 'scalp'})), cfg.tissue = 1; end
if ~isfield(cfg, 'threshold'), cfg.threshold = 0; end
if ~isfield(mri, 'unit'), mri = ft_convert_units(mri); end
  
fprintf('using the segmented MRI\n');

if ~isfield(mri, 'seg') && isequal(cfg.tissue, 1)
  
  if isfield(mri, 'gray') || isfield(mri, 'white') || isfield(mri, 'csf')
    % construct the single image segmentation from the three probabilistic
    % tissue segmentations for csf, white and gray matter
    mri.seg = zeros(size(mri.gray ));
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
      fprintf('smoothing the segmentation with a %d-pixel FWHM kernel\n',cfg.smooth);
      mri.seg = double(mri.seg);
      spm_smooth(mri.seg, mri.seg, cfg.smooth);
    end
    % threshold for the last time
    mri.seg = (mri.seg>(cfg.threshold*max(mri.seg(:))));
  elseif isfield(mri, 'brain')
    mri.seg = mri.brain;
  elseif isfield(mri, 'scalp')
    mri.seg = mri.scalp;
  end
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
  switch mri.unit
    case 'mm'
      scale = scale / 1000;
    case 'cm'
      scale = scale / 100;
    case 'dm'
      scale = scale / 10;
    case 'm'
      scale = scale / 1;
    otherwise
      error('unknown physical dimension in mri.unit');
  end
  if scale~=1
    fprintf('converting MRI surface points from %s into %s\n', cfg.sourceunits, mri.unit);
    pnt = pnt* scale;
  end
  
  bnd(i).pnt = pnt;
  bnd(i).tri = tri;
  fprintf(['segmentation compartment %d of ' num2str(length(cfg.tissue)) ' completed\n'],i);
end

