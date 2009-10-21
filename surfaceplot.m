function [cfg] = surfaceplot(cfg, vol, surf)

% SURFACEPLOT plot functional data on a rendered cortical surface. The
% distance between every voxel in the functional data to each surface point
% on the cortex is computed. Based on this distance, a weighted projection
% of the functional data onto the cortical surface is performed. Activity
% slightly below or otherwise misaligned with the cortical surface is
% therefore also projected onto the surface. The weighted projection also
% smooths the functional data.
%
% Use as
%   surfaceplot(cfg, functional) or
%   surfaceplot(cfg, functional, surface)
% where the functional data is obtained from SOURCEANALYSIS,
% SOURCEINTERPOLATE or VOLUMENORMALIZE.
%
% The configuration can contain
%   cfg.method         = 'nearest' 'sphere_avg', 'sphere_weighteddistance'
%   cfg.sphereradius   = maximum distance from each voxel to the surface to be
%                        included in the sphere projection methods, expressed in mm
%   cfg.surffile       = string, default is 'single_subj_T1.mat' which contains
%                        a triangulation that corresponds with the SPM anatomical
%                        template in MNI coordinates
%   cfg.distmat        = precomputed distance matrix (default = [])
%   cfg.funparameter   = string with the functional parameter of interest
%   cfg.colmin         = functional value mapped to the lowest color (default = 'auto')
%   cfg.colmax         = functional value mapped to the highest color (default = 'auto')
%   cfg.maskparameter  = string with an optional mask parameter
%   cfg.maskcolmin     = mask value mapped to the lowest opacity, i.e. completely transparent (default ='auto')
%   cfg.maskcolmax     = mask value mapped to the highest opacity, i.e. non-transparent (default = 'auto')
%   cfg.downsample     = number (default = 1, i.e. no downsampling)
%   cfg.surfdownsample = number (default = 1, i.e. no downsampling)
%
% The resulting plot can be rotated in 3-D, or you can change the viewpoint
% using the VIEW command.
%
% See also SOURCEPLOT, SLICEINTERP

% Copyright (C) 2006, Jan-Mathijs Schoffelen
% $Log: surfaceplot.m,v $
% Revision 1.7  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.6  2006/07/13 08:48:46  ingnie
% fixed typo's in documentation
%
% Revision 1.5  2006/05/31 10:12:39  roboos
% fixed bug for inside voxels, renamed maxprojdist into sphereradius
%
% Revision 1.4  2006/05/31 07:00:39  roboos
% cleaned up surfaceplot, added options for downsampling volume and surface data, reimplemented nearest neighbour interpolation using distance matrix, flipped dimensions of collin27 surface, added transformation matrix
%
% Revision 1.3  2006/05/29 08:22:28  jansch
% made some changes
%
% Revision 1.2  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.1  2006/01/27 15:41:24  jansch
% first implementation
%

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'funparameter'),   error('cfg.funparameter should be specified'); end
if ~isfield(cfg, 'sphereradius'),   cfg.sphereradius   = []; end
if ~isfield(cfg, 'maskparameter'),  cfg.maskparameter  = []; end
if ~isfield(cfg, 'colmin'),         cfg.colmin         = []; end
if ~isfield(cfg, 'colmax'),         cfg.colmax         = []; end
if ~isfield(cfg, 'maskcolmin'),     cfg.maskcolmin     = []; end
if ~isfield(cfg, 'maskcolmax'),     cfg.maskcolmax     = []; end
if ~isfield(cfg, 'distmat'),        cfg.distmat        = []; end
if ~isfield(cfg, 'downsample'),     cfg.downsample     = 1;  end
if ~isfield(cfg, 'surfdownsample'), cfg.surfdownsample = 1;  end

% downsample the functional data
tmpcfg = [];
tmpcfg.downsample = cfg.downsample;
vol = volumedownsample(tmpcfg, vol);

hasmask = ~isempty(cfg.maskparameter);

if nargin<3
  % default is to use a triangulation that corresponds with the collin27
  % anatomical template in MNI coordinates
  if ~isfield(cfg, 'surffile'),   cfg.surffile = 'single_subj_T1.mat'; end
  tmp = load(cfg.surffile, 'bnd');
  surf = tmp.bnd;
else
  % the cortical surface is supplied by the user
end

if ~isfield(vol, 'transform'),
  vol.transform = eye(4);
end
if ~isfield(surf, 'transform'),
  surf.transform = eye(4);
end

%------extract the stuff that is needed from the input
param      = parameterselection(cfg.funparameter, vol);
dat        = getsubfield(vol, param{1});
if hasmask,
  maskparam  = parameterselection(cfg.maskparameter, vol);
  mask       = getsubfield(vol, maskparam{1});
end
dim        = size(dat);
dimres     = svd(vol.transform(1:3,1:3));
if isfield(vol, 'inside')
  inside = vol.inside;
else
  inside = true(dim);
end
inside = inside(:);

if isempty(cfg.colmin)
  cfg.colmin = min(dat(:));
end
if isempty(cfg.colmax)
  cfg.colmax = max(dat(:));
end

if hasmask && isempty(cfg.maskcolmin)
  cfg.maskcolmin = min(mask(:));
end
if hasmask && isempty(cfg.maskcolmax)
  cfg.maskcolmax = max(mask(:));
end

%------create a matrix, containing the functional voxels in head-coordinates
[X,Y,Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
pos     = [X(:) Y(:) Z(:)]; clear X Y Z;
npos    = size(pos,1);
pos     = warp_apply(vol.transform, pos, 'homogenous');

%------create a matrix, containing the surface positions in head-coordinates
tri  = surf.tri;
pnt  = surf.pnt;
if cfg.surfdownsample>1
  [tri, pnt] = reducepatch(tri, pnt, 1/cfg.surfdownsample);
end
pnt  = warp_apply(surf.transform, pnt, 'homogenous');
npnt = size(pnt,1);

fprintf('%d vertices in surface\n', npnt);
fprintf('%d voxels in functional data\n', npos);

if ~isempty(cfg.distmat)
  %------use the precomputed distance matrix
  distmat = cfg.distmat;
else
  %------compute a distance matrix
  switch cfg.method
    case 'nearest'
      if ~isempty(cfg.sphereradius)
        warning('cfg.sphereradius is not used for method''nearest''');
      end
      % determine the nearest voxel for each surface point
      sub = round(warp_apply(inv(vol.transform), pnt, 'homogenous'));  % express surface vertices in voxel coordinates
      sub(sub(:)<1) = 1;
      sub(sub(:,1)>dim(1),1) = dim(1);
      sub(sub(:,2)>dim(2),2) = dim(2);
      sub(sub(:,3)>dim(3),3) = dim(3);
      ind = sub2ind(dim, sub(:,1), sub(:,2), sub(:,3));
      distmat = sparse(1:npnt, ind, ones(size(ind)), npnt, npos);
      % only voxels inside the brain contain a meaningful functional value
      distmat = distmat(:, inside);
    case {'sphere_avg', 'sphere_weighteddistance', 'sphere_weightedprojection'}
      if isempty(cfg.sphereradius)
        error('cfg.sphereradius should be specified');
      end
      % the distance only has to be computed to voxels inside the brain
      pos  = pos(inside,:);
      npos = size(pos,1);
      % compute the distance between voxels and each surface point
      dpntsq  = sum(pnt.^2,2); % squared distance to origin
      dpossq  = sum(pos.^2,2); % squared distance to origin
      maxnpnt = double(npnt*ceil(4/3*pi*(cfg.sphereradius/max(dimres))^3)); % initial estimate of nonzero entries
      distmat = spalloc(npnt, npos, maxnpnt);
      progress('init', 'textbar', 'computing distance matrix');
      for j = 1:npnt
        progress(j/npnt);
        d   = sqrt(dpntsq(j) + dpossq - 2 * pos * pnt(j,:)');
        sel = find(d<cfg.sphereradius);
        distmat(j, sel) = single(d(sel)) + eps('single');
      end
      progress('close');
    otherwise
      error('unsupported projection method');
  end
end

%------update the configuration
cfg.distmat = distmat;

%------do something with the distance matrix
switch cfg.method
  case 'nearest'
    projmat         = distmat;
  case 'sphere_avg'
    projmat         = distmat;
    [ind1, ind2, d] = find(projmat);
    nnz             = full(sum(spones(projmat),2));
    for k = 1:length(ind1)
      projmat(ind1(k),ind2(k)) = 1./nnz(ind1(k));
    end
  case 'sphere_weighteddistance'
    projmat         = distmat;
    [ind1, ind2, d] = find(projmat);
    projmat         = sparse(ind1, ind2, 1./d, npnt, npos);
    [ind1, ind2, d] = find(projmat);
    normnz          = sqrt(full(sum(projmat.^2, 2)));
    projmat         = sparse(ind1, ind2, d./normnz(ind1), npnt, npos);
  case 'sphere_weightedprojection'
    % JM had something in mind for this, but it is not yet implemented
    error('unsupported projection method');
  otherwise
    error('unsupported projection method');
end

%------compute the values
val = projmat*dat(inside);

if isempty(cfg.colmin), cfg.colmin = min(val(:)); end
if isempty(cfg.colmax), cfg.colmax = max(val(:)); end

%------scale the mask
if hasmask,
  maskval = projmat*mask(inside);
  if isempty(cfg.maskcolmin), cfg.maskcolmin = min(maskval(:)); end
  if isempty(cfg.maskcolmax), cfg.maskcolmax = max(maskval(:)); end
  %if ~isempty(cfg.maskcolmax), maskval(find(maskval>cfg.maskcolmax)) = cfg.maskcolmax; end
  %if ~isempty(cfg.maskcolmin), maskval(find(maskval<cfg.maskcolmin)) = cfg.maskcolmin; end
  %maskval = (maskval-min(maskval(:)))./max(maskval(:));
  %scaling can be done by setting the alim property.
end

%------do the plotting
brain = repmat([0.781 0.762 0.664], length(val), 1);
h1 = patch('Vertices', pnt, 'Faces', tri, 'FaceVertexCData', brain , 'FaceColor', 'interp');
set(h1, 'EdgeColor', 'none');
axis   off;
axis vis3d;
axis equal;

h2 = patch('Vertices', pnt, 'Faces', tri, 'FaceVertexCData', val , 'FaceColor', 'interp');
set(h2, 'EdgeColor', 'none');
if hasmask
  set(h2, 'FaceVertexAlphaData', maskval);
  set(h2, 'FaceAlpha',          'interp');
  set(h2, 'AlphaDataMapping',   'scaled');
  alim(gca, [cfg.maskcolmin cfg.maskcolmax]);
end
caxis(gca,[cfg.colmin cfg.colmax]);

lighting gouraud
colormap jet
camlight

