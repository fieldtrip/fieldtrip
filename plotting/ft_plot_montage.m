function ft_plot_montage(dat, varargin)

transform = keyval('transform',   varargin);
loc       = keyval('location',    varargin);
ori       = keyval('orientation', varargin);
mask      = keyval('datmask',     varargin);
cmap      = keyval('colormap',    varargin);
srange    = keyval('slicerange',  varargin);
nslice    = keyval('nslice',      varargin);

% set the location if empty
if isempty(loc) && (isempty(transform) || all(all(transform-eye(4)==0)==1))
  % go to the middle of the volume if the data seem to be in voxel coordinates
  loc = size(dat)./2;
elseif isempty(loc)
  % otherwise take the origin of the coordinate system
  loc = [0 0 0];
end

% check compatibility of inputs
if size(loc, 1) == 1 && isempty(nslice)
  nslice = 20;
elseif size(loc, 1) == 1 && ~isempty(nslice)
  % this is not a problem, slice spacing will be determined
elseif size(loc, 1) > 1 && isempty(nslice)
  % this is not a problem, number of slices is determined by loc
  nslice = size(loc, 1);
elseif size(loc, 1) > 1 && ~isempty(nslice)
  if size(loc, 1) ~= nslice
    error('you should either specify a set of locations or a single location with a number of slices');
  end
end

% set the transformation matrix if empty
if isempty(transform)
  transform = eye(4);
end

% set the orientation if empty
if isempty(ori),
  ori = [0 0 1];
end

% ensure the ori to have unit norm
for k = 1:size(ori,1)
  ori(k,:) = ori(k,:)./norm(ori(k,:));
end

% determine the slice range
if size(loc, 1) == 1 && nslice > 1,
  if isempty(srange) || (ischar(srange) && strcmp(srange, 'auto'))
    srange = [-50 70];
  else
  end
  loc = repmat(loc, [nslice 1]) + linspace(srange(1),srange(2),nslice)'*ori;    
end

% ensure that the ori has the same size as the loc
if size(ori,1)==1 && size(loc,1)>1,
  ori = repmat(ori, size(loc,1), 1);
end

div    = [ceil(sqrt(nslice)) ceil(sqrt(nslice))];
for k = 1:nslice
  ix     = mod(k-1, div(1));
  iy     = floor((k-1)/div(1));
  h(k)   = ft_plot_slice(dat, 'transform', transform, 'location', loc(k,:), ...
                              'orientation', ori(k,:), 'colormap', cmap, 'datmask', mask);
  siz    = size(get(h(k), 'xdata'));
  set(h(k), 'xdata', ix*(siz(1) -1) + repmat((0:siz(1)-1)', [1 siz(2)]));
  set(h(k), 'ydata', iy*(siz(2) -1) + repmat((0:siz(2)-1) , [siz(1) 1]));
  set(h(k), 'zdata', zeros(siz));   
end
axis equal;
axis tight;
axis off;