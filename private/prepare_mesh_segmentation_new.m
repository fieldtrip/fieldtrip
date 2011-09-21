function bnd = prepare_mesh_segmentation(mri,varargin)

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

% process the inputs
dim = mri.dim;
tissuelabels = ft_getopt('tissuelabels',varargin,[]);
threshold    = ft_getopt('threshold',varargin,0);
smooth       = ft_getopt('smooth',varargin,'no');
numvertices  = ft_getopt('smooth',varargin,[]);
tissuetype   = ft_getopt('tissuetype',varargin,[]);
unit         = ft_getopt('unit',varargin,'mm');
if ~isfield(mri, 'unit'), mri = ft_convert_units(mri); end
% tissueval    = ft_getopt('smooth',varargin,[]);% not sure about it


  if isfield(mri, 'gray') || isfield(mri, 'white') || isfield(mri, 'csf')
    % construct the single image segmentation from the three probabilistic
    % tissue segmentations for csf, white and gray matter
    mri.seg = zeros(size(mri.gray )); % FIXME
    if isfield(mri, 'gray')
      fprintf('including gray matter in segmentation for brain compartment\n')
      mri.seg = mri.seg | (mri.gray>(threshold*max(mri.gray(:)))); % FIXME
    end
    if isfield(mri, 'white')
      fprintf('including white matter in segmentation for brain compartment\n')
      mri.seg = mri.seg | (mri.white>(threshold*max(mri.white(:)))); % FIXME
    end
    if isfield(mri, 'csf')
      fprintf('including CSF in segmentation for brain compartment\n')
      mri.seg = mri.seg | (mri.csf>(threshold*max(mri.csf(:)))); % FIXME
    end
    if ~strcmp(cfg.smooth, 'no'),
      fprintf('smoothing the segmentation with a %d-pixel FWHM kernel\n',cfg.smooth);
      mri.seg = double(mri.seg);
      spm_smooth(mri.seg, mri.seg, smooth); % FIXME
    end
    % threshold for the last time
    mri.seg = (mri.seg>(cfg.threshold*max(mri.seg(:)))); % FIXME
  elseif isfield(mri, 'brain')
    mri.seg = mri.brain; % FIXME
  elseif isfield(mri, 'scalp')
    mri.seg = mri.scalp; % FIXME
  end

  [mrix, mriy, mriz] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));

  % construct the triangulations of the boundaries from the segmented MRI
  for i=1:length(cfg.tissuetype) % FIXME
    fprintf('triangulating the boundary of compartment %d\n', i);
    seg = imfill((mri.seg==tissueval(i)), 'holes'); % FIXME
    ori(1) = mean(mrix(seg(:)));
    ori(2) = mean(mriy(seg(:)));
    ori(3) = mean(mriz(seg(:)));
    [pnt, tri] = triangulate_seg(seg, numvertices(i), ori); % is tri okay?
    %     tri = projecttri(pnt);
    % apply the coordinate transformation from voxel to head coordinates
    pnt = warp_apply(mri.transform,pnt);
    % rescale
    pnt = scaleunit(unit,pnt);
    % output
    bnd(i).pnt = pnt;
    bnd(i).tri = tri;
    fprintf(['segmentation compartment %d of ' num2str(length(cfg.tissue)) ' completed\n'],i);
  end

function pnt = scaleunit(unit,pnt)
  % convert the MRI surface points into the same units as the source/gradiometer
  scale = 1;
  switch unit
    case 'mm'
      scale = scale * 1000;
    case 'cm'
      scale = scale * 100;
    case 'dm'
      scale = scale * 10;
    case 'm'
      scale = scale * 1;
    otherwise
      error('unknown physical dimension in input ''unit''');
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
