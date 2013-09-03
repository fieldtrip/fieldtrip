function mesh=prepare_mesh_hexahedral(cfg,mri)

% PREPARE_MESH_HEXAHEDRAL
%
% See also PREPARE_MESH_SEGMENTATION, PREPARE_MESH_MANUAL, PREPARE_MESH_HEADSHAPE
%
% Configuration options for generating a regular 3-D grid
%   cfg.tissue = cell with the names of the compartments that should be
%   meshed
%   cfg.resolution = desired resolution of the mesh (standard = 1)
%
% Copyrights (C) 2012-2013, Johannes Vorwerk
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% get the default options
cfg.tissue      = ft_getopt(cfg, 'tissue');
cfg.resolution  = ft_getopt(cfg, 'resolution');
cfg.shift       = ft_getopt(cfg, 'shift');
cfg.background  = ft_getopt(cfg, 'background');

if isempty(cfg.tissue)
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  fn = fieldnames(mri);
  for i=1:numel(fn),if numel(mri.(fn{i}))==prod(mri.dim), segfield=fn{i};end;end
  cfg.tissue=setdiff(unique(mri.(segfield)(:)),0);
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
else
  % the code below assumes that it is an indexed representation
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
end

if isempty(cfg.resolution)
    warning('Using standard resolution 1 mm')
    cfg.resolution = 1;
end

if isempty(cfg.shift)
    warning('No node-shift selected')
    cfg.shift = 0;
elseif cfg.shift > 0.3
    warning('Node-shift should not be larger than 0.3')
    cfg.shift = 0.3;
end

if isempty(cfg.background)
    cfg.background = 0;
end

% do the mesh extraction
% this has to be adjusted for FEM!!!
if iscell(cfg.tissue)
  % this assumes that it is a probabilistic representation
  % for example {'brain', 'skull', scalp'}
  try
    temp = zeros(size(mri.(cfg.tissue{1})(:)));
    for i=1:numel(cfg.tissue)
      temp = [temp,mri.(cfg.tissue{i})(:)];
    end
    [val,seg] = max(temp,[],2);
    seg = seg - 1;
    seg = reshape(seg,mri.dim);
  catch
    error('Please specify cfg.tissue to correspond to tissue types in the segmented MRI')
  end
  tissue = cfg.tissue;
else
  % this assumes that it is an indexed representation
  % for example [3 2 1]
  seg = zeros(mri.dim);
  tissue = {};
  for i=1:numel(cfg.tissue)
    seg = seg + i*(mri.seg==cfg.tissue(i));
    if isfield(mri, 'seglabel')
      try
        tissue{i} = mri.seglabel{cfg.tissue(i)};
      catch
        error('Please specify cfg.tissue to correspond to (the name or number of) tissue types in the segmented MRI')
      end
    else
      tissue{i} = sprintf('tissue %d', i);
    end
  end
end

% reslice to desired resolution

if (cfg.resolution ~= 1)
    % this should be done like this: split seg into probabilistic, reslice
    % single compartments, take maximum values
    seg_array = [];
    
    seg_indices = unique(seg);    
    
    for i=1:(length(unique(seg)))
        seg_reslice.anatomy = double(seg == (i-1));
        seg_reslice.dim = mri.dim;
        seg_reslice.transform = eye(4);
        seg_reslice.transform(1:3,4) = -ceil(mri.dim/2);
        
        cfg_reslice = [];
        cfg_reslice.resolution = cfg.resolution;
        cfg_reslice.dim = ceil(mri.dim/cfg.resolution);
        
        seg_build = ft_volumereslice(cfg_reslice,seg_reslice);
        
        seg_array = [seg_array,seg_build.anatomy(:)];
        
        clear seg_reslice;
    end
    
    [max_seg seg_build.seg] = max(seg_array,[],2);
    
    clear max_seg seg_array;
    
    seg_build.seg = reshape(seg_build.seg,seg_build.dim);
    seg_build.seg = seg_indices(seg_build.seg);
    seg_build.transform = mri.transform;
    
    clear seg_build.anatomy;
else
    seg_build.seg = seg;
    seg_build.dim = mri.dim;
    
    clear seg;
end

% ensure that the segmentation is binary and that there is a single contiguous region
% FIXME is this still needed when it is already binary?
%seg = volumethreshold(seg, 0.5, tissue);

ft_hastoolbox('simbio', 1);

% build the mesh

mesh = build_mesh_hexahedral(seg_build,cfg);

% converting position of meshpoints to the head coordinate system

if (cfg.resolution ~= 1)
    mesh.pnt = cfg.resolution * mesh.pnt;
end

mesh.pnt = warp_apply(mri.transform,mesh.pnt,'homogeneous');

labels = mesh.labels;

clear mesh.labels;

mesh.tissue = zeros(size(labels));
numlabels = size(unique(labels),1);
mesh.tissuelabel = {};
ulabel = sort(unique(labels));
for i = 1:numlabels
  mesh.tissue(labels == ulabel(i)) = i;
  mesh.tissuelabel{i} = tissue{i};
end

end % function
