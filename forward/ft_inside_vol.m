function [inside] = ft_inside_vol(dippos, headmodel, varargin)

% FT_INSIDE_VOL locates dipole locations inside/outside the source
% compartment of a volume conductor model.
%
% Use as
%   [inside] = ft_inside_vol(dippos, headmodel, ...)
%
% The input should be
%   dippos      = Nx3 matrix with dipole positions
%   headmodel   = structure with volume conductor model
% and the output is
%   inside      = boolean vector indicating for each dipole wether it is inside the source compartment
%
% Additional optional input arguments should be given in key value pairs and can include
%   inwardshift = number
%   grad        = structure with gradiometer information, used for localspheres
%   headshape   = structure with headshape, used for old CTF localspheres strategy

% Copyright (C) 2003-2013, Robert Oostenveld
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

% get the optional input arguments
grad        = ft_getopt(varargin, 'grad');
headshape   = ft_getopt(varargin, 'headshape');
inwardshift = ft_getopt(varargin, 'inwardshift');

% for backward compatibility
headmodel = fixpos(headmodel);

% determine the type of volume conduction model
switch ft_voltype(headmodel)
  
  case {'singlesphere' 'concentricspheres'}
    if ~isfield(headmodel, 'source')
      % locate the innermost compartment and remember it
      [dum, headmodel.source] = min(headmodel.r);
    end
    if isfield(headmodel, 'o')
      % shift dipole positions toward origin of sphere
      tmp = dippos - repmat(headmodel.o, size(dippos,1), 1);
    else
      tmp = dippos;
    end
    distance = sqrt(sum(tmp.^2, 2))-headmodel.r(headmodel.source);
    % positive if outside, negative if inside
    inside   = distance<0;
    
  case 'localspheres'
    if ~isempty(headshape) && ~isempty(grad)
      % use the specified headshape to construct the bounding triangulation
      [pos, tri] = headsurface(headmodel, grad, 'headshape', headshape, 'inwardshift', inwardshift, 'surface', 'skin');
      inside = bounding_mesh(dippos, pos, tri);
    elseif ~isempty(grad)
      % use the volume conductor model to construct an approximate headshape
      [pos, tri] = headsurface(headmodel, grad, 'inwardshift', inwardshift, 'surface', 'skin');
      inside = bounding_mesh(dippos, pos, tri);
    else
      % only check whether the dipole is in any of the spheres
      nspheres = size(headmodel.r,1);
      ndipoles = size(dippos,1);
      inside = zeros(ndipoles,1);
      for sph=1:nspheres
        % temporary shift dipole positions toward origin
        if isfield(headmodel, 'o')
          tmp = dippos - repmat(headmodel.o(sph,:), [ndipoles 1]);
        else
          tmp = dippos;
        end
        flag = (sqrt(sum(tmp.^2,2)) <= headmodel.r(sph));
        inside = inside + flag;
      end
      inside = inside>0;
    end
    
  case {'infinite' 'infinite_monopole'}
    % an empty headmodel in combination with gradiometers indicates a magnetic dipole
    % in an infinite vacuum, i.e. all dipoles can be considered to be inside
    inside = true(1,size(dippos,1));
    
  case {'halfspace', 'halfspace_monopole'}
    inside = false(1,size(dippos,1));
    for i = 1:size(dippos,1);
      pol = dippos(i,:);
      % condition of dipoles/monopoles falling in the non conductive halfspace
      inside(i) = acos(dot(headmodel.ori,(pol-headmodel.pos)./norm(pol-headmodel.pos))) >= pi/2;
    end
    
  case 'slab_monopole'
    inside = false(1,size(dippos,1));
    for i=1:size(dippos,1);
      pol = dippos(i,:);
      % condition of dipoles/monopoles falling in the non conductive halfspace
      % Attention: voxels on the boundary are automatically considered outside the strip
      instrip1  = acos(dot(headmodel.ori1,(pol-headmodel.pos1)./norm(pol-headmodel.pos1))) > pi/2;
      instrip2  = acos(dot(headmodel.ori2,(pol-headmodel.pos2)./norm(pol-headmodel.pos2))) > pi/2;
      inside(i) = instrip1 & instrip2;
    end
    
  case {'bem', 'dipoli', 'bemcp', 'openmeeg', 'asa', 'singleshell', 'neuromag'}
    % this is a model with a realistic shape described by a triangulated boundary
    [pos, tri] = headsurface(headmodel, [], 'inwardshift', inwardshift, 'surface', 'brain');
    inside = bounding_mesh(dippos, pos, tri);
    
  case {'simbio'}
    
    % FIXME we have to rethink which tissue types should be flagged as inside
    tissue = intersect({'gray', 'white', 'csf', 'brain'}, headmodel.tissuelabel);
    brain = false(size(headmodel.tissue));
    for i=1:numel(tissue)
      fprintf('selecting source positions inside "%s"\n', tissue{i});
      brain = brain | (headmodel.tissue == find(strcmp(headmodel.tissuelabel, tissue{i})));
    end
     
    % prune the mesh, i.e. only retain hexaheders that are part of the brain
    discard = ~brain;
    
    fprintf('pruning mesh from %d to %d elements (%d%%)\n', length(discard), sum(discard), round(100*sum(~discard)/length(discard)));
    
    subhex = headmodel.hex(~discard,:);
    
    % determine the center of each volume element
    elementpos = zeros(size(subhex,1),3);
    for i=1:8
      elementpos = elementpos + headmodel.pos(subhex(:,i),:);
    end
    elementpos = elementpos/8;
    
    % find the nearest hexaheder for each of the source points
    % FIXME this is only guaranteed to work for regular hexahedral meshes
    stopwatch = tic;
    subindx = dsearchn(elementpos, dippos(1,:));
    t = toc(stopwatch);
    fprintf('determining inside points, this takes about %d seconds\n', round(size(dippos,1)*t));
    subindx = dsearchn(elementpos, dippos);
    
    % The following code is only guaranteed to work with convex elements. Regular
    % hexahedra and tetrahedra are convex, and the adapted hexahedra we can use with
    % SIMBIO must also be convex.
    inside = false(size(subindx));
    for i=1:numel(subindx)
      % determine for each source position whether it is actually inside the nearest hexaheder
      posx = headmodel.pos(subhex(subindx(i),:),1);
      posy = headmodel.pos(subhex(subindx(i),:),2);
      posz = headmodel.pos(subhex(subindx(i),:),3);
      pos = [posx posy posz];
      inside(i) = isequal(convhull(pos), convhull([pos; dippos(i,:)]));
    end
        
  otherwise
    error('unrecognized volume conductor model');
end

% ensure that it is a boolean column vector
inside(isnan(inside(:))) = 0;
inside = logical(inside(:));
