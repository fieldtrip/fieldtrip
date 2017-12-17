function varargout = interp_gridded(transform, val, pos_to, varargin)

% INTERP_GRIDDED computes a matrix that interpolates values that were
% observed on positions in a regular 3-D grid onto positions that are
% unstructured, e.g. the vertices of a cortical sheet.
%
% Use as
%   [val]                = interp_gridded(transform, val, pos, ...) or
%   [interpmat, distmat] = interp_gridded(transform, val, pos, ...)
% where
%   transform  homogenous coordinate transformation matrix for the volume
%   val        3-D matrix with the values in the volume
%   pos        Mx3 matrix with the vertex positions onto which the data should
%              be interpolated
% 
% Optional arguments are specified in key-value pairs and can be
%    projmethod   = 'nearest', 'sphere_avg', 'sphere_weighteddistance'
%    sphereradius = number
%    distmat      = NxM matrix with precomputed distances
%    inside       = indices for inside voxels (or logical array)

% Copyright (C) 2007-2015, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
  ft_error('Not enough input arguments.');
end

% get the optional arguments
projvec      = ft_getopt(varargin, 'projvec',       1);
projweight   = ft_getopt(varargin, 'projweight',    ones(size(projvec)));
projcomb     = ft_getopt(varargin, 'projcomb',      'mean'); % or max
projthresh   = ft_getopt(varargin, 'projthresh',    []);
projmethod   = ft_getopt(varargin, 'projmethod');    % required
sphereradius = ft_getopt(varargin, 'sphereradius');  % required for some projection methods
distmat      = ft_getopt(varargin, 'distmat');       % will be computed if not present
inside       = ft_getopt(varargin, 'inside');

dim    = size(val);
dimres = svd(transform(1:3,1:3)); % to reduce the number of elements in the distance matrix
npos_to   = size(pos_to,1);
npos_from = prod(dim);

if isempty(distmat)
  % compute the distance matrix
  switch projmethod
    case 'nearest'
      % determine the nearest voxel for each vertex
      sub = round(ft_warp_apply(inv(transform), pos_to, 'homogenous'));  % express
      sub(sub(:)<1) = 1;
      sub(sub(:,1)>dim(1),1) = dim(1);
      sub(sub(:,2)>dim(2),2) = dim(2);
      sub(sub(:,3)>dim(3),3) = dim(3);
      ind = sub2ind(dim, sub(:,1), sub(:,2), sub(:,3));
      distmat = sparse(1:npos_to, ind, ones(size(ind)), npos_to, npos_from);
      if ~isempty(inside)
        % only voxels inside the brain contain a meaningful functional value
        distmat = distmat(:, inside);
      end

    case {'sphere_avg', 'sphere_weighteddistance', 'sphere_weightedprojection'}
      if isempty(sphereradius)
        ft_error('sphereradius should be specified');
      end

      [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
      pos_from  = ft_warp_apply(sparse(transform), [X(:) Y(:) Z(:)]);
      % the distance only has to be computed to voxels inside the brain
      pos_from  = pos_from(inside,:);
      npos_from = size(pos_from,1);
      % compute the distance between all voxels and each surface point
      dfromsq = sum(pos_from.^2,2); % squared distance to origin
      dtosq   = sum(pos_to.^2,  2); % squared distance to origin
      maxnpnt = double(npos_to*ceil(4/3*pi*(sphereradius/max(dimres))^3)); % initial estimate of nonzero entries
      distmat = spalloc(npos_to, npos_from, maxnpnt);
      ft_progress('init', 'textbar', 'computing distance matrix');
      for j = 1:npos_to
        ft_progress(j/npos_to);
        d   = sqrt(dfromsq + dtosq(j) - 2 * pos_from * pos_to(j,:)');
        sel = find(d<sphereradius);
        distmat(j, sel) = single(d(sel)) + eps('single');
      end
      ft_progress('close');

    case 'project'
      % do nothing I believe
    
    otherwise
      ft_error('unsupported projection method');
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
    projmat         = sparse(ind1, ind2, 1./d, npos_to, npos_from);
    [ind1, ind2, d] = find(projmat);
    normnz          = sqrt(full(sum(projmat.^2, 2)));
    projmat         = sparse(ind1, ind2, d./normnz(ind1), npos_to, npos_from);
  
  case 'project'
      % this method is Joachim's implementation that was originally in
      % ft_sourceplot, it assumes the functional data to be defined on a
      % regular 3D grid, and that the transformation to world-space is known
      
      % we also need the dim
      dim = ft_getopt(varargin, 'dim');
      dat = zeros(size(pos_to,1),1);
      
      % convert projvec in mm to a factor, assume mean distance of 70mm
      projvec = (70-projvec)/70;
      for iproj = 1:length(projvec),
        sub = round(ft_warp_apply(inv(transform), pos_to*projvec(iproj), 'homogenous'));  % express
        sub(sub(:)<1) = 1;
        sub(sub(:,1)>dim(1),1) = dim(1);
        sub(sub(:,2)>dim(2),2) = dim(2);
        sub(sub(:,3)>dim(3),3) = dim(3);
        ind = sub2ind(dim, sub(:,1), sub(:,2), sub(:,3));
        if strcmp(projcomb,'mean')
          dat = dat + projweight(iproj) * val(ind);
        elseif strcmp(projcomb,'max')
          dat  = max([dat projweight(iproj) * val(ind)],[],2);
          tmp2 = min([dat projweight(iproj) * val(ind)],[],2);
          fi   = find(dat < max(tmp2));
          val(fi) = tmp2(fi);
        else
          ft_error('undefined method to combine projections; use cfg.projcomb= mean or max')
        end
      end
      if strcmp(projcomb,'mean'),
        dat = dat/length(projvec);
      end
      %     if ~isempty(projthresh),
      %       mm=max(abs(val(:)));
      %       maskval(abs(val) < projthresh*mm) = 0;
      %     end
      %
    
  otherwise
    ft_error('unsupported projection method');
end  % case projmethod

if nargout==1 && ~strcmp(projmethod, 'project')
  % return the interpolated values
  varargout{1} = projmat * val(:);
elseif nargout==1
  varargout{1} = dat;
else
  % return the interpolation and the distance matrix
  varargout{1} = projmat;
  varargout{2} = distmat;
end
