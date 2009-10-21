function varargout = interp_gridded(transform, val, pnt, varargin)

% INTERP_GRIDDED computes a matrix that interpolates values that were
% observed on a regular 3-D grid onto a random set of points.
%
% Use as
%   [val]                = interp_gridded(transform, val, pnt, ...) or
%   [interpmat, distmat] = interp_gridded(transform, val, pnt, ...)
% where
%   transform  homogenous coordinate transformation matrix for the volume
%   val        3-D matrix with the values in the volume
%   pnt        Mx3 matrix with the vertex positions onto which the data should
%              be interpolated
% 
% Optional arguments are specified in key-value pairs and can be
%    projmethod   = 'nearest', 'sphere_avg', 'sphere_weighteddistance'
%    sphereradius = number
%    distmat      = NxM matrix with precomputed distances
%    inside       = indices for inside voxels (or logical array)

% Copyright (C) 2007, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: interp_gridded.m,v $
% Revision 1.1  2007/05/06 12:07:26  roboos
% new function, based on old surfaceplot and interp_ungridded
% this implementation is much faster for sourceplot with method=surface
%
% Revision 1.1  2007/02/07 07:42:06  roboos
% new implementation to be used in teh new sourceplot, mainly based on code from the old surfaceplot
%

if nargin<3
  error('Not enough input arguments.');
end

% get the optional arguments
projmethod   = keyval('projmethod',    varargin);   % required
sphereradius = keyval('sphereradius',  varargin);   % required for some projection methods
distmat      = keyval('distmat',       varargin);   % will be computed if not present
inside       = keyval('inside',        varargin);

dim = size(val);
dimres = svd(transform(1:3,1:3)); % to reduce the number of elements in the distance matrix
npnt = size(pnt,1);
npos = prod(dim);

if isempty(distmat)
  % compute the distance matrix
  switch projmethod
    case 'nearest'
      % determine the nearest voxel for each vertex
      sub = round(warp_apply(inv(transform), pnt, 'homogenous'));  % express
      sub(sub(:)<1) = 1;
      sub(sub(:,1)>dim(1),1) = dim(1);
      sub(sub(:,2)>dim(2),2) = dim(2);
      sub(sub(:,3)>dim(3),3) = dim(3);
      ind = sub2ind(dim, sub(:,1), sub(:,2), sub(:,3));
      distmat = sparse(1:npnt, ind, ones(size(ind)), npnt, npos);
      if ~isempty(inside)
        % only voxels inside the brain contain a meaningful functional value
        distmat = distmat(:, inside);
      end

    case {'sphere_avg', 'sphere_weighteddistance', 'sphere_weightedprojection'}
      if isempty(sphereradius)
        error('sphereradius should be specified');
      end

      [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
      pos = warp_apply(sparse(transform), [X(:) Y(:) Z(:)]);
      % the distance only has to be computed to voxels inside the brain
      pos  = pos(inside,:);
      npos = size(pos,1);
      % compute the distance between all voxels and each surface point
      dpntsq  = sum(pnt.^2,2); % squared distance to origin
      dpossq  = sum(pos.^2,2); % squared distance to origin
      maxnpnt = double(npnt*ceil(4/3*pi*(sphereradius/max(dimres))^3)); % initial estimate of nonzero entries
      distmat = spalloc(npnt, npos, maxnpnt);
      progress('init', 'textbar', 'computing distance matrix');
      for j = 1:npnt
        progress(j/npnt);
        d   = sqrt(dpntsq(j) + dpossq - 2 * pos * pnt(j,:)');
        sel = find(d<sphereradius);
        distmat(j, sel) = single(d(sel)) + eps('single');
      end
      progress('close');

    otherwise
      error('unsupported projection method');
  end % case projmethod
end % if isempty distmat

%------do something with the distance matrix
switch projmethod
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
    projmat         = sparse(ind1, ind2, 1./d, npnt, npnt1);
    [ind1, ind2, d] = find(projmat);
    normnz          = sqrt(full(sum(projmat.^2, 2)));
    projmat         = sparse(ind1, ind2, d./normnz(ind1), npnt, npnt1);

  otherwise
    error('unsupported projection method');
end  % case projmethod

if nargout==1
  % return the interpolated values
  varargout{1} = projmat * val(:);
else
  % return the interpolation and the distance matrix
  varargout{1} = projmat;
  varargout{2} = distmat;
end

