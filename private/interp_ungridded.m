function varargout = interp_ungridded(pntin, pntout, varargin)

% INTERP_UNGRIDDED computes an interpolation matrix for two
% clouds of 3-D points
%
% Use as
%   [val] = interp_ungridded(pntin, pntout, 'data', valin, ...)
% where
%   pntin   Nx3 matrix with the vertex positions
%   pntout  Mx3 matrix with the vertex positions onto which the data should
%           be interpolated
%   valin   NxK matrix with functional data, defined at the points in pntin
%
% Alternatively to get the interpolation matrix itself, you can use it as
%   [interpmat, distmat] = interp_ungridded(pntin, pntout, ...)
% 
%
% Optional arguments are specified in key-value pairs and can be
%    projmethod   = 'nearest', 'sphere_avg', 'sphere_weighteddistance',
%                   'smudge'
%    sphereradius = number
%    distmat      = NxM matrix with precomputed distances
%    triin        = triangulation for the first set of vertices
%    triout       = triangulation for the second set of vertices
%    data         = NxK matrix with functional data

% Copyright (C) 2007-2011, Jan-Mathijs Schoffelen & Robert Oostenveld
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
projmethod   = ft_getopt(varargin, 'projmethod');     % required
sphereradius = ft_getopt(varargin, 'sphereradius');   % required for some projection methods
distmat      = ft_getopt(varargin, 'distmat');        % will be computed if needed and not present
triin        = ft_getopt(varargin, 'triin');          % not yet implemented
triout       = ft_getopt(varargin, 'triout');   
val          = ft_getopt(varargin, 'data');           % functional data defined at pntin

hasval  = ~isempty(val);
npntin  = size(pntin, 1);
npntout = size(pntout, 1);

dimres  = sqrt(sum((pntout(2,:)-pntout(1,:)).^2,2));

if isempty(distmat)
  %------compute a distance matrix
  switch projmethod
    case 'nearest'
      if ~isempty(sphereradius)
        warning('sphereradius is not used for projmethod ''nearest''');
      end
      % determine the nearest voxel for each surface point
      ind     = find_nearest(pntout, pntin, 5);
      distmat = sparse(1:npntout, ind, ones(size(ind)), npntout, npntin);

    case {'sphere_avg', 'sphere_weighteddistance'}
      if isempty(sphereradius)
        error('sphereradius should be specified');
      end
      % compute the distance between voxels and each surface point
      dpntinsq  = sum(pntin.^2,2); % squared distance to origin
      dpntoutsq  = sum(pntout.^2,2); % squared distance to origin
      maxnpnt = double(npntout*ceil(4/3*pi*(sphereradius/max(dimres))^3)); % initial estimate of nonzero entries
      maxnpnt = min(maxnpnt, npntout*npntin);
      %ft_progress('init', 'none', 'computing distance matrix');
      val = nan(maxnpnt, 1);
      indx1 = nan(maxnpnt, 1);
      indx2 = nan(maxnpnt, 1);
      cnt = 1;
      for j = 1:npntout
        %ft_progress(j/npntout);
        %d   = dpntoutsq(j) + dpntinsq - 2 * pntin * pntout(j,:)';
        %sel = find(d<sphereradius.^2);
        
        % the following lines are equivalent to the previous 2 but use
        % fewer flops
        d    = sqrt(dpntinsq - 2 * pntin * pntout(j,:)' + dpntoutsq(j));
        sel  = find(d < sphereradius);
        
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
      distmat = sparse(indx1, indx2, val, npntout, npntin);
      %ft_progress('close');

    case 'smudge'
      if isempty(triout),
        error('the ''smudge'' method needs a triangle definition');
      end
      [datin, loc] = ismember(pntout, pntin, 'rows');
      [datout, S1] = smudge(datin, triout, 6); %FIXME 6 is number of iterations, improve here
    
      sel = find(datin);
      S2  = sparse(sel(:), loc(datin), ones(npntin,1), npntout, npntin);
      distmat = S1 * S2;
    
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
    projmat         = sparse(ind1, ind2, 1./d, npnt, npntin);
    [ind1, ind2, d] = find(projmat);
    normnz          = sqrt(full(sum(projmat.^2, 2)));
    projmat         = sparse(ind1, ind2, d./normnz(ind1), npnt, npntin);

  case 'smudge'
    projmat = distmat;
    
  otherwise
    error('unsupported projection method');
end  % case projmethod

if hasval
  % return the interpolated values
  varargout{1} = projmat * val;
else
  % return the interpolation and the distance matrix
  varargout{1} = projmat;
  varargout{2} = distmat;
end

