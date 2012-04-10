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
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<3
  error('Not enough input arguments.');
end

% get the optional arguments
projmethod   = ft_getopt(varargin, 'projmethod');    % required
sphereradius = ft_getopt(varargin, 'sphereradius');  % required for some projection methods
distmat      = ft_getopt(varargin, 'distmat');       % will be computed if not present
inside       = ft_getopt(varargin, 'inside');

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
      ft_progress('init', 'textbar', 'computing distance matrix');
      for j = 1:npnt
        ft_progress(j/npnt);
        d   = sqrt(dpntsq(j) + dpossq - 2 * pos * pnt(j,:)');
        sel = find(d<sphereradius);
        distmat(j, sel) = single(d(sel)) + eps('single');
      end
      ft_progress('close');

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

