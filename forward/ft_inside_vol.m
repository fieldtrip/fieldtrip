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
    
  case {'simbio'}
    
    brain = false(size(vol.tissue));
    brain = brain | vol.tissue == find(strcmp(vol.tissuelabel, 'gray'));
    brain = brain | vol.tissue == find(strcmp(vol.tissuelabel, 'white'));
    brain = brain | vol.tissue == find(strcmp(vol.tissuelabel, 'csf'));
    
    minbrain = min(vol.pos(vol.hex(brain(:)), :), [], 1);
    maxbrain = max(vol.pos(vol.hex(brain(:)), :), [], 1);
    
    minpos = min(pos, [], 1);
    maxpos = max(pos, [], 1);
    
    % combine the two bounding boxes
    minbox = max([minbrain; minpos], [], 1);
    maxbox = min([maxbrain; maxpos], [], 1);
    
    % prune the mesh to the bounding box
    discard1 = true(size(vol.hex,1),1);
    discard2 = true(size(vol.hex,1),1);
    for i=1:8
      discard1 = discard1 & any(bsxfun(@minus, vol.pos(vol.hex(:,i),:), minbox)<0,2);
      discard2 = discard2 & any(bsxfun(@minus, vol.pos(vol.hex(:,i),:), maxbox)>0,2);
    end
    discard = discard1 | discard2;
    
    fprintf('pruning mesh from %d to %d elements (%d%%)\n', length(discard), sum(discard), round(100*sum(discard)/length(discard)));
    
    vol.hex    = vol.hex(~discard,:);
    vol.tissue = vol.tissue(~discard);
    
    % determine the center of each volume element
    elementpos = zeros(size(vol.hex,1),3);
    for i=1:8
      elementpos = elementpos + vol.pos(vol.hex(:,i),:);
    end
    elementpos = elementpos/8;
    
    stopwatch = tic;
    k = dsearchn(elementpos, pos(1,:));
    t = toc(stopwatch);
    fprintf('determining inside points, this takes about %d seconds\n', round(size(pos,1)*t));
    k = dsearchn(elementpos, pos);
    
    % select the source positions that are inside the brain
    inside = brain(k);
    
  otherwise
    error('unrecognized volume conductor model');
end

% ensure that it is a boolean column vector
inside(isnan(inside(:))) = 0;
inside = logical(inside(:));
