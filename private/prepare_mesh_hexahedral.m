function mesh = prepare_mesh_hexahedral(cfg, mri)

% PREPARE_MESH_HEXAHEDRAL
%
% Configuration options for generating a regular 3-D grid
%   cfg.tissue     = cell with the names of the compartments that should be meshed
%   cfg.resolution = desired resolution of the mesh (default = 1)
%   cfg.shift
%   cfg.background
%
% See also PREPARE_MESH_SEGMENTATION, PREPARE_MESH_MANUAL, PREPARE_MESH_HEADSHAPE

% Copyrights (C) 2012-2013, Johannes Vorwerk
%
% $Id$

% ensure that the input is consistent with what this function expects
mri = ft_checkdata(mri, 'datatype', {'volume', 'segmentation'}, 'hasunit', 'yes');

% get the default options
cfg.tissue      = ft_getopt(cfg, 'tissue');
cfg.resolution  = ft_getopt(cfg, 'resolution', 1);  % this is in mm
cfg.shift       = ft_getopt(cfg, 'shift');
cfg.background  = ft_getopt(cfg, 'background', 0);

if isempty(cfg.tissue)
  mri = ft_datatype_segmentation(mri, 'segmentationstyle', 'indexed');
  fn = fieldnames(mri);
  for i = 1:numel(fn),
    if (numel(mri.(fn{i})) == prod(mri.dim)) && (~strcmp(fn{i}, 'inside'))
      segfield = fn{i};
    end
  end
  cfg.tissue = setdiff(unique(mri.(segfield)(:)), 0);
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

if isempty(cfg.shift)
  warning('No node-shift selected')
  cfg.shift = 0;
elseif cfg.shift > 0.3
  warning('Node-shift should not be larger than 0.3')
  cfg.shift = 0.3;
end


% do the mesh extraction
% this has to be adjusted for FEM!!!
if iscell(cfg.tissue)
  % this assumes that it is a probabilistic representation
  % for example {'brain', 'skull', scalp'}
  try
    temp = zeros(size(mri.(cfg.tissue{1})(:)));
    for i = 1:numel(cfg.tissue)
      temp = [temp, mri.(cfg.tissue{i})(:)];
    end
    [val, seg] = max(temp, [], 2);
    seg = seg - 1;
    seg = reshape(seg, mri.dim);
  catch
    error('Please specify cfg.tissue to correspond to tissue types in the segmented MRI')
  end
  tissue = cfg.tissue;
else
  % this assumes that it is an indexed representation
  % for example [3 2 1]
  seg = zeros(mri.dim);
  tissue = {};
  for i = 1:numel(cfg.tissue)
    seg = seg + i*(mri.seg == cfg.tissue(i));
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
  % this should be done like this: split seg into probabilistic, reslice single compartments, take maximum values
  seg_array = [];
  
  seg_indices = unique(seg);
  
  for i = 1:(length(unique(seg)))
    seg_reslice.anatomy   = double(seg == (i-1));
    seg_reslice.dim       = mri.dim;
    seg_reslice.transform = eye(4);
    seg_reslice.transform(1:3, 4) = -ceil(mri.dim/2);
    
    cfg_reslice = [];
    cfg_reslice.resolution = cfg.resolution;
    cfg_reslice.dim = ceil(mri.dim/cfg.resolution);
    
    seg_build = ft_volumereslice(cfg_reslice, seg_reslice);
    
    seg_array = [seg_array, seg_build.anatomy(:)];
    
    clear seg_reslice;
  end
  
  [max_seg seg_build.seg] = max(seg_array, [], 2);
  
  clear max_seg seg_array;
  
  seg_build.seg = reshape(seg_build.seg, seg_build.dim);
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

mesh = build_mesh_hexahedral(cfg, seg_build);

% converting position of meshpoints to the head coordinate system

if (cfg.resolution ~= 1)
  mesh.pnt = cfg.resolution * mesh.pnt;
end

mesh.pnt = ft_warp_apply(mri.transform, mesh.pnt, 'homogeneous');

labels = mesh.labels;

clear mesh.labels;

mesh.tissue = zeros(size(labels));
numlabels = size(unique(labels), 1);
mesh.tissuelabel = {};
ulabel = sort(unique(labels));
for i = 1:numlabels
  mesh.tissue(labels == ulabel(i)) = i;
  mesh.tissuelabel{i} = tissue{i};
end


end % function

%% subfunctions %%

function mesh = build_mesh_hexahedral(cfg, mri)

background = cfg.background;
shift = cfg.shift;
% extract number of voxels in each direction
% x_dim = mri.dim(1);
% y_dim = mri.dim(2);
% z_dim = mri.dim(3);
%labels = mri.seg;
fprintf('Dimensions of the segmentation before restriction to bounding-box: %i %i %i\n', mri.dim(1), mri.dim(2), mri.dim(3));

[bb_x, bb_y, bb_z] = ind2sub(size(mri.seg), find(mri.seg));
shift_coord = [min(bb_x) - 2, min(bb_y) - 2, min(bb_z) - 2];
bb_x = [min(bb_x), max(bb_x)];
bb_y = [min(bb_y), max(bb_y)];
bb_z = [min(bb_z), max(bb_z)];
x_dim = size(bb_x(1)-1:bb_x(2)+1, 2);
y_dim = size(bb_y(1)-1:bb_y(2)+1, 2);
z_dim = size(bb_z(1)-1:bb_z(2)+1, 2);
labels = zeros(x_dim, y_dim, z_dim);
labels(2:(x_dim-1), 2:(y_dim-1), 2:(z_dim-1)) = mri.seg(bb_x(1):bb_x(2), bb_y(1):bb_y(2), bb_z(1):bb_z(2));

fprintf('Dimensions of the segmentation after restriction to bounding-box: %i %i %i\n', x_dim, y_dim, z_dim);

% create elements

mesh.hex = create_elements(x_dim, y_dim, z_dim);
fprintf('Created elements...\n' )


% create nodes

mesh.pnt = create_nodes(x_dim, y_dim, z_dim);
fprintf('Created nodes...\n' )


if(shift < 0 | shift > 0.3)
  error('Please choose a shift parameter between 0 and 0.3!');
elseif(shift > 0)
  
  mesh.pnt = shift_nodes(mesh.pnt, mesh.hex, labels, shift, x_dim, y_dim, z_dim);
  
end

%background = 1;
% delete background voxels(if desired)
if(background == 0)
  mesh.hex = mesh.hex(labels ~= 0, :);
  mesh.labels = labels(labels ~= 0);
else
  mesh.labels = labels(:);
end


% delete unused nodes
[C, ia, ic] = unique(mesh.hex(:));
mesh.pnt = mesh.pnt(C, :, :, :);
mesh.pnt = mesh.pnt + repmat(shift_coord, size(mesh.pnt, 1), 1);
mesh.hex(:) = ic;

end % subfunction

% function creating elements from a MRI-Image with the dimensions x_dim,
% y_dim, z_dim. Each voxel of the MRI-Image corresponds to one element in
% the hexahedral mesh. The numbering of the elements is as follows:
% the first x-y-plane(z == 1) is numbered by incrementing in x-direction
% first until x_dim is reached. This is done for each row in y-direction
% until y_dim is reached. Using the resulting x_dim-by-y_dim numbering as
% an overlay and adding (i-1)*(x_dim*y_dim) to the overlay while i is
% iterating through the z-dimension(until i == z_dim) we obtain a numbering
% for all the elements.
% The node-numbering is done accordingly: the bottom left node in element i
% has number i in the node-numbering. All the remaining nodes are numbered
% in the same manner as described above for the element numbering(note the
% different dimensionalities: x_dim+1 instead of x_dim etc.).
function elements = create_elements(x_dim, y_dim, z_dim)
elements = zeros(x_dim*y_dim*z_dim, 8);
% create an offset vector for the bottom-left nodes in each element
b = 1:((x_dim+1)*(y_dim));
% delete the entries where the node does not correspond to an element's
% bottom-left node(i.e. where the x-component of the node is equal to
% x_dim+1)
b = b(mod(b, (x_dim+1)) ~= 0);
% repeat offset to make it fit the number of elements
b = repmat(b, 1, (z_dim));

% create vector accounting for the offset of the nodes in z-direction
c = fix([0:((x_dim)*(y_dim)*(z_dim)-1)]/((x_dim)*(y_dim))) * (x_dim+1)*(y_dim+1);

% create the elements by assigning the nodes to them. entries 1 through 4
% describe the bottom square of the hexahedron, entries 5 through 8 the top
% square.
elements(:, 1) = b + c;
elements(:, 2) = b + c + 1;
elements(:, 3) = b + c + (x_dim+1) + 1;
elements(:, 4) = b + c + (x_dim+1);
elements(:, 5) = b + c + (x_dim + 1)*(y_dim+1);
elements(:, 6) = b + c + (x_dim + 1)*(y_dim+1) + 1;
elements(:, 7) = b + c + (x_dim + 1)*(y_dim+1) + (x_dim+1) + 1;
elements(:, 8) = b + c + (x_dim + 1)*(y_dim+1) + (x_dim+1);
clear b;
clear c;
end %subfunction

% function creating the nodes and assigning coordinates in the
% [0, x_dim]x[0, y_dim]x[0, z_dim] box. for details on the node-numbering see
% comments for create_elements.
function nodes = create_nodes(x_dim, y_dim, z_dim)
nodes = zeros(((x_dim+1)*(y_dim+1)*(z_dim + 1)), 3);
% offset vector for node coordinates
b = [0:((x_dim+1)*(y_dim+1)*(z_dim+1)-1)];

% assign coordinates within the box
nodes(:, 1) = mod(b, (x_dim+1));
nodes(:, 2) = mod(fix(b/(x_dim+1)), (y_dim+1));
nodes(:, 3) = fix(b/((x_dim + 1)*(y_dim+1)));

clear b;

end

% function shifting the nodes
function nodes = shift_nodes(points, hex, labels, sh, x_dim, y_dim, z_dim)
cfg = [];
fprintf('Applying shift %f\n', sh);
nodes = points;

% helper vector for indexing the nodes
b = 1:(x_dim+1)*(y_dim+1)*(z_dim+1);

% vector which gives the corresponding element to a node(in the sense
% of how the node-numbering was done, see comment for create_elements).
% note that this vector still contains values for nodes which do not have a
% corresponding element. this will be manually omitted in the finding
% of the surrounding elements(see below).
offset = (b - nodes(b, 2)' - (nodes(b, 3))'*(y_dim+1+x_dim))';
offset(offset <= 0) = size(labels, 1)+1;
offset(offset > size(hex, 1)) = size(labels, 1)+1;

% create array containing the surrounding elements for each node
%surrounding = zeros((x_dim+1)*(y_dim+1)*(z_dim+1), 8);

% find out the surrounding of each node(if there is any)
% find element to which the node is bottom left front
%surrounding((nodes(b, 1) < x_dim) & (nodes(b, 3) < z_dim), 1) = offset((nodes(b, 1) < x_dim) & (nodes(b, 3) < z_dim));
% bottom right front
%surrounding(nodes(b, 1) > 0, 2) = offset(nodes(b, 1) > 0, 1) -1;
% bottom left back
%surrounding(nodes(b, 2) > 0, 3) = offset(nodes(b, 2) > 0, 1) - x_dim;
% bottom right back
%surrounding((nodes(b, 2) > 0) & (nodes(b, 1) > 0), 4) = offset((nodes(b, 2) > 0) & (nodes(b, 1) > 0), 1) - (x_dim) - 1;
% top left front
%surrounding(nodes(b, 3) > 0, 5) = offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim);
% top right front
%surrounding(nodes(b, 3) > 0, 6) = offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim) - 1;
% top left back
%surrounding(nodes(b, 3) > 0, 7) = offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim) - x_dim;
% top right back
%surrounding((nodes(b, 3) > 0) & (nodes(b, 2) > 0), 8) = offset((nodes(b, 3) > 0) & (nodes(b, 2) > 0), 1) - (x_dim)*(y_dim) - (x_dim) -1;
%clear offset;
%clear b;

% those entries in the surrounding matrix which point to non-existing
% elements(> size(hex, 1) or <= 0) we overwrite with a
% dummy element
%surrounding(surrounding <= 0) = size(labels, 1) + 1;
%surrounding(surrounding > size(hex, 1)) = size(labels, 1)+1;

% set the label of the dummy element to be zero(background)
labels(size(labels, 1)+1) = 0;

% matrixs holding the label of each surrounding element
%surroundinglabels = labels(surrounding);
surroundinglabels = zeros((x_dim+1)*(y_dim+1)*(z_dim+1), 8);
surroundinglabels((nodes(b, 1) < x_dim) & (nodes(b, 3) < z_dim), 1) = labels(offset((nodes(b, 1) < x_dim) & (nodes(b, 3) < z_dim)));
surroundinglabels(nodes(b, 1) > 0, 2) = labels(offset(nodes(b, 1) > 0, 1) -1);
surroundinglabels(nodes(b, 2) > 0, 3) = labels(offset(nodes(b, 2) > 0, 1) - x_dim);
offsetnow = offset((nodes(b, 2) > 0) & (nodes(b, 1) > 0), 1) - (x_dim) - 1;
offsetnow(offsetnow <= 0) = size(labels, 1)+1;
surroundinglabels((nodes(b, 2) > 0) & (nodes(b, 1) > 0), 4) = labels(offsetnow);
offsetnow = offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim);
offsetnow(offsetnow <= 0) = size(labels, 1)+1;
surroundinglabels(nodes(b, 3) > 0, 5) = labels(offsetnow);
offsetnow = offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim) - 1;
offsetnow(offsetnow <= 0) = size(labels, 1)+1;
surroundinglabels(nodes(b, 3) > 0, 6) = labels(offsetnow);
offsetnow = offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim) - x_dim;
offsetnow(offsetnow <= 0) = size(labels, 1)+1;
surroundinglabels(nodes(b, 3) > 0, 7) = labels(offsetnow);
offsetnow = offset((nodes(b, 3) > 0) & (nodes(b, 2) > 0), 1) - (x_dim)*(y_dim) - (x_dim) -1;
offsetnow(offsetnow <= 0) = size(labels, 1)+1;
surroundinglabels((nodes(b, 3) > 0) & (nodes(b, 2) > 0), 8) = labels(offsetnow);

% matrix showing how many types of each label are around a given node
distribution = zeros(size(nodes, 1), size(unique(labels), 1));
% the following assumes that the labels are 1, 2, ...
for l = 1:size(unique(labels), 1)
  for k = 1:8
    distribution(:, l) = distribution(:, l) + (surroundinglabels(:, k) == l);
  end
end

% fill up the last column with the amount of background labels
distribution(:, (size(unique(labels), 1))) = 8 - sum(distribution(:, 1:size(unique(labels), 1))');

% how many different labels are there around each node
distsum = sum(distribution>0, 2);

% set zeros to Inf in order to make the finding of a minimum
% meaningful
distribution(distribution == 0) = Inf;

% find out which is the minority label
[mins, minpos] = min(distribution, [], 2);
clear distribution;

% calculate the centroid for each element
centroids = zeros(size(hex, 1), 3);
for l = 1:3
  centroids(:, l) = sum(reshape(nodes(hex(:, :), l), size(hex, 1), 8)')'/8;
end

% set a dummy centroid
centroids(size(centroids, 1) +1, :) = 0;

tbc = zeros(size(surroundinglabels));
% helper matrix, c(i, j, k) is one when surroundinglabels(i, j) == k
for i = 1:size(unique(labels), 1)+1
  c = zeros(size(surroundinglabels, 2), size(unique(labels), 1)+1);
  if (i == size(unique(labels), 1)+1)
    c = surroundinglabels == 0;
  else
    c = surroundinglabels == i;
  end
  tbc(ismember(minpos, i) == 1, :) = c(ismember(minpos, i) == 1, :);
end

%     % matrix that shows which elements are to be considered as minority
%     % around a given node
%     clear surroundinglabels;
%      % helper matrix, c(i, j, k) is one when surroundinglabels(i, j) == k
%     c = zeros(size(surroundinglabels, 1), size(surroundinglabels, 2), size(unique(labels), 1)+1);
%     for i = 1:size(unique(labels), 1)
%         c(:, :, i) = surroundinglabels == i;
%     end
%     c(:, :, size(unique(labels), 1)+1) = surroundinglabels == 0;
%
%
%     % matrix that shows which elements are to be considered as minority
%     % around a given node
%     tbc = zeros(size(surroundinglabels));
%     clear surroundinglabels;
%     for i = 1:size(unique(labels), 1)+1
%     tbc(ismember(minpos, i) == 1, :) = c(ismember(minpos, i) == 1, :, i);
%     end
clear c;

% delete cases in which we don't have a real minimum
tbcsum = sum(tbc, 2);
tbc(tbcsum == 8, :) = 0;
tbc(tbcsum == 4, :) = 0;
tbcsum((distsum>2) & (mins > 1), :) = 0;
tbcsum(tbcsum == 8) = 0;
tbcsum(tbcsum == 4) = 0;
tbcsum((distsum>2) & (mins > 1)) = 0;

%surroundingconsidered = surrounding.*tbc;
surroundingconsidered = zeros(size(tbc));
surroundingconsidered((nodes(b, 1) < x_dim) & (nodes(b, 3) < z_dim), 1) = offset((nodes(b, 1) < x_dim) & (nodes(b, 3) < z_dim)).*tbc((nodes(b, 1) < x_dim) & (nodes(b, 3) < z_dim), 1);
surroundingconsidered(nodes(b, 1) > 0, 2) = (offset(nodes(b, 1) > 0, 1) -1).*(tbc(nodes(b, 1) > 0, 2));
surroundingconsidered(nodes(b, 2) > 0, 3) = (offset(nodes(b, 2) > 0, 1) - x_dim).*tbc(nodes(b, 2) > 0, 3);
surroundingconsidered((nodes(b, 2) > 0) & (nodes(b, 1) > 0), 4) = (offset((nodes(b, 2) > 0) & (nodes(b, 1) > 0), 1) - (x_dim) - 1).*tbc((nodes(b, 2) > 0) & (nodes(b, 1) > 0), 4).*tbc((nodes(b, 2) > 0) & (nodes(b, 1) > 0), 4);
surroundingconsidered(nodes(b, 3) > 0, 5) = (offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim)).*tbc(nodes(b, 3) > 0, 5);
surroundingconsidered(nodes(b, 3) > 0, 6) = (offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim) - 1).*tbc(nodes(b, 3) > 0, 6);
surroundingconsidered(nodes(b, 3) > 0, 7) = (offset(nodes(b, 3) > 0, 1) - (x_dim)*(y_dim) - x_dim).*tbc(nodes(b, 3) > 0, 7);
surroundingconsidered((nodes(b, 3) > 0) & (nodes(b, 2) > 0), 8) = (offset((nodes(b, 3) > 0) & (nodes(b, 2) > 0), 1) - (x_dim)*(y_dim) - (x_dim) -1).*tbc((nodes(b, 3) > 0) & (nodes(b, 2) > 0), 8);
%clear surrounding;
clear tbc;

tbcsum(tbcsum == 8) = 0;
tbcsum(tbcsum == 4) = 0;
tbcsum((distsum>2) & (mins > 1)) = 0;

clear distsum;
% get the surrounding elements which are to be considered for the shift


% use the dummy centroid to make computations easier
surroundingconsidered(surroundingconsidered == 0) = size(centroids, 1);

% calculate the combination of the centroids which are to be considered
% for the shift
centroidcomb = zeros(size(nodes));
centroidcomb(:, 1) = sum(reshape(centroids(surroundingconsidered, 1), [], 8), 2);
centroidcomb(:, 2) = sum(reshape(centroids(surroundingconsidered, 2), [], 8), 2);
centroidcomb(:, 3) = sum(reshape(centroids(surroundingconsidered, 3), [], 8), 2);
clear surroundingconsidered;
centroidcomb(tbcsum ~= 0, 1) = centroidcomb(tbcsum ~= 0, 1)./tbcsum(tbcsum ~= 0);
centroidcomb(tbcsum ~= 0, 2) = centroidcomb(tbcsum ~= 0, 2)./tbcsum(tbcsum ~= 0);
centroidcomb(tbcsum ~= 0, 3) = centroidcomb(tbcsum ~= 0, 3)./tbcsum(tbcsum ~= 0);

% finally apply the shift
nodes(tbcsum == 0, :) = points(tbcsum == 0, :);
nodes(tbcsum ~= 0, :) = (1-sh)*nodes(tbcsum ~= 0, :) + sh*centroidcomb(tbcsum ~= 0, :);

end %subfunction



