function ft_plot_montage(dat, varargin)

% FT_PLOT_MONTAGE

loc       = ft_getopt(varargin, 'location');
ori       = ft_getopt(varargin, 'orientation');
srange    = ft_getopt(varargin, 'slicerange');
slicesize = ft_getopt(varargin, 'slicesize');
nslice    = ft_getopt(varargin, 'nslice');
transform = ft_getopt(varargin, 'transform', eye(4));

dointersect = ~isempty(ft_getopt(varargin, 'intersectmesh'));

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
optarg = varargin;
for k = 1:nslice
  
  % define 'x' and 'y' axis in projection plane, the definition of x and y is more or less arbitrary
  [x, y] = projplane(ori(k,:)); % z = ori
  
  % get the transformation matrix to project onto the xy-plane
  T  = [x(:) y(:) ori(k,:)' loc(k,:)'; 0 0 0 1];
  
  optarg = ft_setopt(optarg, 'location',    loc(k,:));
  optarg = ft_setopt(optarg, 'orientation', ori(k,:));
  ix     = mod(k-1, div(1));
  iy     = floor((k-1)/div(1));
  h(k)   = ft_plot_slice(dat, optarg{:});
  
  xtmp = get(h(k), 'xdata');
  ytmp = get(h(k), 'ydata');
  ztmp = get(h(k), 'zdata');
  siz  = size(xtmp);
  if k==1 && isempty(slicesize)
    slicesize = siz;
  end
  
  % project the positions onto the xy-plane
  pos = [xtmp(:) ytmp(:) ztmp(:)];
  pos = warp_apply(inv(T), pos);
  
  xtmp = reshape(pos(:,1), siz);
  ytmp = reshape(pos(:,2), siz);
  ztmp = reshape(pos(:,3), siz);
  
  % add some offset in the x and y directions to create the montage
  offset(1) = ix*(slicesize(1)-1);
  offset(2) = iy*(slicesize(2)-1); 
  
  % update the positions
  set(h(k), 'xdata', offset(1) + ytmp);
  set(h(k), 'ydata', offset(2) + xtmp);
  set(h(k), 'zdata',         0 * ztmp);
  
  if dointersect,
    if ~exist('pprevious', 'var'), pprevious = []; end
    p = setdiff(findobj(gcf, 'type', 'patch'), pprevious);
    for kk = 1:numel(p)
      xtmp = get(p(kk), 'xdata');
      ytmp = get(p(kk), 'ydata');
      ztmp = get(p(kk), 'zdata');
      siz2 = size(xtmp);
      
      pos = [xtmp(:) ytmp(:) ztmp(:)];
      pos = warp_apply(inv(T), pos);
  
      xtmp = reshape(pos(:,1), siz2);
      ytmp = reshape(pos(:,2), siz2);
      ztmp = reshape(pos(:,3), siz2);
  
      % update the positions
      set(p(kk), 'xdata', offset(1) + ytmp);
      set(p(kk), 'ydata', offset(2) + xtmp);
      set(p(kk), 'zdata',         0 * ztmp);
    end
    pprevious = [pprevious(:);p(:)];
  end
  drawnow;
end
set(gcf, 'color', [0 0 0]);
set(gca, 'zlim', [0 1]);
%axis equal;
axis off;
view([0 90]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y] = projplane(z)
[u, s, v] = svd([eye(3) z(:)]);
x = u(:, 2)';
y = u(:, 3)';
