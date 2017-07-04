function varargout = interp_ungridded(pos_from, pos_to, varargin)

% INTERP_UNGRIDDED computes an interpolation matrix for two
% clouds of 3-D points
%
% Use as
%   [val] = interp_ungridded(pos_from, pos_to, 'data', valin, ...)
% where
%   pos_from  Nx3 matrix with the vertex positions
%   pos_to    Mx3 matrix with the vertex positions onto which the data should
%             be interpolated
%
% Alternatively to get the interpolation matrix itself, you can use it as
%   [interpmat, distmat] = interp_ungridded(pos_from, pos_to, ...)
% 
%
% Optional arguments are specified in key-value pairs and can be
%    projmethod   = 'nearest', 'sphere_avg', 'sphere_weighteddistance',
%                   'smudge'
%    sphereradius = number
%    distmat      = NxM matrix with precomputed distances
%    triout       = triangulation for the second set of vertices
%    data         = NxK matrix with functional data
%
% Optional extra arguments when using projmethod = 'sphere_weighteddistance'
%    power        = power parameter as in the Inverse Distance Weighting
%                   function proposed by Shepard (default = 1).

% Copyright (C) 2007-2013, Jan-Mathijs Schoffelen & Robert Oostenveld
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
projmethod   = ft_getopt(varargin, 'projmethod');     % required
sphereradius = ft_getopt(varargin, 'sphereradius');   % required for some projection methods
powerparam   = ft_getopt(varargin, 'power', 1);
distmat      = ft_getopt(varargin, 'distmat');        % will be computed if needed and not present
triout       = ft_getopt(varargin, 'triout');   
dat          = ft_getopt(varargin, 'data');           % functional data defined at pos_from
inside       = ft_getopt(varargin, 'inside');

hasdat    = ~isempty(dat);
hasinside = ~isempty(inside);
npos_from = size(pos_from, 1);
npos_to   = size(pos_to, 1);

dimres  = sqrt(sum((pos_to(2,:)-pos_to(1,:)).^2,2));

if hasinside,
  % convert to boolean vector
  tmp         = false(npos_from,1);
  tmp(inside) = true;
  inside      = tmp;
  clear tmp;
else
  inside      = true(npos_from,1);
end

if isempty(distmat)
  %------compute a distance matrix
  switch projmethod
    case 'nearest'
      if ~isempty(sphereradius)
        ft_warning('sphereradius is not used for projmethod ''nearest''');
      end
      % determine the nearest voxel for each surface point
      ind     = find_nearest(pos_to, pos_from, 5);
      sel     = ind>0;
      indx    = 1:npos_to;
      distmat = sparse(indx(sel), ind(sel), ones(size(ind(sel))), npos_to, npos_from);

    case {'sphere_avg', 'sphere_weighteddistance'}
      if isempty(sphereradius)
        ft_error('sphereradius should be specified');
      end
      % compute the distance between voxels and each surface point
      dpos_fromsq = sum(pos_from.^2,2); % squared distance to origin
      dpos_tosq   = sum(pos_to.^2,2); % squared distance to origin
      maxnpnt = double(npos_to*ceil(4/3*pi*(sphereradius/max(dimres))^3)); % initial estimate of nonzero entries
      maxnpnt = min(maxnpnt, npos_to*npos_from);
      %ft_progress('init', 'none', 'computing distance matrix');
      val   = nan(maxnpnt, 1);
      indx1 = nan(maxnpnt, 1);
      indx2 = nan(maxnpnt, 1);
      cnt = 1;
      for j = 1:npos_to
        %ft_progress(j/npos_to);
        %d   = dpos_tosq(j) + dpos_fromsq - 2 * pos_from * pos_to(j,:)';
        %sel = find(d<sphereradius.^2);
        
        % the following lines are equivalent to the previous 2 but use
        % fewer flops
        d    = sqrt(dpos_fromsq - 2 * pos_from * pos_to(j,:)' + dpos_tosq(j));
        sel  = find(d < sphereradius & inside);
        
        nsel = numel(sel);
        if nsel>0
        indx1(cnt:(cnt+nsel-1)) = j(ones(nsel,1));
        indx2(cnt:(cnt+nsel-1)) = sel(:);
        val(cnt:(cnt+nsel-1))   = d(sel) + eps('double');
        cnt = cnt + nsel;
        end
      end
      indx1(isnan(indx1)) = [];
      indx2(isnan(indx2)) = [];
      val(isnan(val))     = [];
      distmat = sparse(indx1, indx2, val, npos_to, npos_from);
      %ft_progress('close');

    case 'smudge'
      if isempty(triout),
        ft_error('the ''smudge'' method needs a triangle definition');
      end
      [datin, loc] = ismember(pos_to, pos_from, 'rows');
      [datout, S1] = smudge(datin, triout, 6); %FIXME 6 is number of iterations, improve here
    
      sel = find(datin);
      S2  = sparse(sel(:), loc(datin), ones(npos_from,1), npos_to, npos_from);
      distmat = S1 * S2;
    
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
    projmat         = sparse(ind1, ind2, d.^-powerparam, npos_to, npos_from);
    [ind1, ind2, d] = find(projmat);
    normnz          = full(sum(projmat, 2));
    projmat         = sparse(ind1, ind2, d./normnz(ind1), npos_to, npos_from);

  case 'smudge'
    projmat = distmat;
  
 
  otherwise
    ft_error('unsupported projection method');
end  % case projmethod

if hasdat
  % return the interpolated values
  varargout{1} = projmat * dat;
else
  % return the interpolation and the distance matrix
  varargout{1} = projmat;
  varargout{2} = distmat;
end

