function [inside] = ft_inside_vol(pos, vol)

% FT_INSIDE_VOL locates dipole locations inside/outside the source
% compartment of a volume conductor model.
%
% [inside] = ft_inside_vol(pos, vol, ...)
%
% where the input should be
%   pos      Nx3 matrix with dipole positions
%   vol      structure with volume conductor model
% and the output is
%   inside   list of dipoles inside the brain compartment
%            (1=inside, 0=outisde)
%
% Additional optional input arguments should be given in key value pairs
% and can include
%   <none>

% Copyright (C) 2003-2007, Robert Oostenveld
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

% determine the type of volume conduction model
switch ft_voltype(vol)

  % single-sphere or multiple concentric spheres
  case {'singlesphere' 'concentricspheres'}
    if ~isfield(vol, 'source')
      % locate the innermost compartment and remember it
      [dum, vol.source] = min(vol.r);
    end
    if isfield(vol, 'o')
      % shift dipole positions toward origin of sphere
      tmp = pos - repmat(vol.o, size(pos,1), 1);
    else
      tmp = pos;
    end
    distance = sqrt(sum(tmp.^2, 2))-vol.r(vol.source);
    % positive if outside, negative if inside
    inside   = distance<0;

    % multiple overlapping sphere volume conductor model
  case 'localspheres'

    % nspheres = size(vol.r,1);
    % ndipoles = size(pos,1);
    % inside = zeros(ndipoles,1);
    % for sph=1:nspheres
    % for dip=1:ndipoles
    %   if inside(dip)
    %     % the dipole has already been detected in one of the other spheres
    %     continue
    %   end
    %   inside(dip) = (norm(pos(dip,:) - vol.o(sph,:)) <= vol.r(sph));
    % end
    % end
    % outside = find(inside==0);
    % inside  = find(inside==1);

    % this is a much faster implementation
    nspheres = size(vol.r,1);
    ndipoles = size(pos,1);
    inside = zeros(ndipoles,1);
    for sph=1:nspheres
      % temporary shift dipole positions toward origin
      if isfield(vol, 'o')
        tmp = pos - repmat(vol.o(sph,:), [ndipoles 1]);
      else
        tmp = pos;
      end
      flag = (sqrt(sum(tmp.^2,2)) <= vol.r(sph));
      inside = inside + flag;
    end
    inside  = inside>0;

    % model with a realistic shape described by a triangulated boundary
  case {'bem', 'dipoli', 'bemcp', 'asa', 'singleshell', 'neuromag'}
    if ~isfield(vol, 'source')
      % locate the innermost compartment and remember it
      vol.source = find_innermost_boundary(vol.bnd);
    end
    % use the specified source compartment
    pnt = vol.bnd(vol.source).pnt;
    tri = vol.bnd(vol.source).tri;
    % determine the dipole positions that are inside the brain compartment
    inside  = bounding_mesh(pos, pnt, tri);

    % unrecognized volume conductor model
  otherwise
    error('unrecognized volume conductor model');
end

% ensure that these are column vectors
inside(find(isnan(inside(:)))) = 0;
inside = logical(inside(:));
