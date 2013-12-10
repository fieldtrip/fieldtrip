function [inside] = ft_inside_vol(pos, vol, varargin)

% FT_INSIDE_VOL locates dipole locations inside/outside the source
% compartment of a volume conductor model.
%
% Use as
%   [inside] = ft_inside_vol(pos, vol, ...)
%
% The input should be
%   pos         = Nx3 matrix with dipole positions
%   vol         = structure with volume conductor model
% and the output is
%   inside      = boolean vector indicating for each dipole wether it is inside the source compartment
%
% Additional optional input arguments should be given in key value pairs and can include
%   inwardshift = number
%   grad        = structure with gradiometer information, used for localspheres
%   headshape   = structure with headshape, used for old CTF localspheres strategy

% Copyright (C) 2003-2013, Robert Oostenveld
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

% get the optional input arguments
grad        = ft_getopt(varargin, 'grad');
headshape   = ft_getopt(varargin, 'headshape');
inwardshift = ft_getopt(varargin, 'inwardshift');

% determine the type of volume conduction model
switch ft_voltype(vol)
  
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
    
  case 'localspheres'
    if ~isempty(headshape) && ~isempty(grad)
      % use the specified headshape to construct the bounding triangulation
      [pnt, tri] = headsurface(vol, grad, 'headshape', headshape, 'inwardshift', inwardshift, 'surface', 'skin');
      inside = bounding_mesh(pos, pnt, tri);
    elseif ~isempty(grad)
      % use the volume conductor model to construct an approximate headshape
      [pnt, tri] = headsurface(vol, grad, 'inwardshift', inwardshift, 'surface', 'skin');
      inside = bounding_mesh(pos, pnt, tri);
    else
      % only check whether the dipole is in any of the spheres
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
      inside = inside>0;
    end
    
  case {'infinite' 'infinite_monopole'}
    % an empty vol in combination with gradiometers indicates a magnetic dipole
    % in an infinite vacuum, i.e. all dipoles can be considered to be inside
    inside = true(1,size(pos,1));
    
  case {'halfspace', 'halfspace_monopole'}
    inside = false(1,size(pos,1));
    for i = 1:size(pos,1);
      pol = pos(i,:);
      % condition of dipoles/monopoles falling in the non conductive halfspace
      inside(i) = acos(dot(vol.ori,(pol-vol.pnt)./norm(pol-vol.pnt))) >= pi/2;
    end
    
  case 'slab_monopole'
    inside = false(1,size(pos,1));
    for i=1:size(pos,1);
      pol = pos(i,:);
      % condition of dipoles/monopoles falling in the non conductive halfspace
      % Attention: voxels on the boundary are automatically considered outside the strip
      instrip1  = acos(dot(vol.ori1,(pol-vol.pnt1)./norm(pol-vol.pnt1))) > pi/2;
      instrip2  = acos(dot(vol.ori2,(pol-vol.pnt2)./norm(pol-vol.pnt2))) > pi/2;
      inside(i) = instrip1 & instrip2;
    end
    
  case {'bem', 'dipoli', 'bemcp', 'openmeeg', 'asa', 'singleshell', 'neuromag'}
    % this is a model with a realistic shape described by a triangulated boundary
    [pnt, tri] = headsurface(vol, [], 'inwardshift', inwardshift, 'surface', 'brain');
    inside = bounding_mesh(pos, pnt, tri);
    
  otherwise
    error('unrecognized volume conductor model');
end

% ensure that it is a boolean column vector
inside(isnan(inside(:))) = 0;
inside = logical(inside(:));
