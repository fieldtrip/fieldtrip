function [inside] = ft_inside_headmodel(dippos, headmodel, varargin)

% FT_INSIDE_HEADMODEL locates dipole locations inside/outside the source
% compartment of a volume conductor model.
%
% Use as
%   [inside] = ft_inside_headmodel(dippos, headmodel, ...)
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
%   surface     = string, the compartment that dipoles should be inside, either
%                 'brain' (default) or 'scalp'/'skin'. Choose 'scalp'/'skin' to
%                 fit muscle and eye sources outside the brain, as in the HArtMuT model.

% Copyright (C) 2003-2024, Robert Oostenveld
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
surface     = ft_getopt(varargin, 'surface', 'brain');

% for backward compatibility
headmodel = fixpos(headmodel);

% determine the type of volume conduction model
switch ft_headmodeltype(headmodel)

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
      inside = surface_inside(dippos, pos, tri);
    elseif ~isempty(grad)
      % use the volume conductor model to construct an approximate headshape
      [pos, tri] = headsurface(headmodel, grad, 'inwardshift', inwardshift, 'surface', 'skin');
      inside = surface_inside(dippos, pos, tri);
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
    for i = 1:size(dippos,1)
      pol = dippos(i,:);
      % condition of dipoles/monopoles falling in the non conductive halfspace
      inside(i) = acos(dot(headmodel.ori,(pol-headmodel.pos)./norm(pol-headmodel.pos))) >= pi/2;
    end

  case 'slab_monopole'
    inside = false(1,size(dippos,1));
    for i=1:size(dippos,1)
      pol = dippos(i,:);
      % condition of dipoles/monopoles falling in the non conductive halfspace
      % Attention: voxels on the boundary are automatically considered outside the strip
      instrip1  = acos(dot(headmodel.ori1,(pol-headmodel.pos1)./norm(pol-headmodel.pos1))) > pi/2;
      instrip2  = acos(dot(headmodel.ori2,(pol-headmodel.pos2)./norm(pol-headmodel.pos2))) > pi/2;
      inside(i) = instrip1 & instrip2;
    end

  case {'bem', 'dipoli', 'bemcp', 'openmeeg', 'asa', 'singleshell', 'neuromag', 'nolte', 'hbf'}
    % this is a model with a realistic shape described by a triangulated boundary
    switch surface
      case 'brain'
        [pos, tri] = headsurface(headmodel, [], 'inwardshift', inwardshift, 'surface', 'brain');
        % surface_inside returns 1, 0, or nan for a degenerate point that lands on the mesh,
        % the ==1 coerces it to a clean logical so the caller always gets a boolean vector
        inside = surface_inside(dippos, pos, tri) == 1;
      case {'scalp', 'skin'}
        % a dipole is in the scalp compartment when it is inside the skin but outside the skull
        if ~isfield(headmodel, 'bnd') || numel(headmodel.bnd)<2
          ft_error('the scalp compartment requires a head model with at least two boundaries');
        end
        [skinpos, skintri] = headsurface(headmodel, [], 'inwardshift', inwardshift, 'surface', 'skin');
        order  = surface_nesting(headmodel.bnd, 'outsidefirst');
        skull  = order(2); % the boundary just inside the skin
        % coerce each test to a logical first, a degenerate point gives nan and nan cannot
        % be combined with the & operator
        inskin  = surface_inside(dippos, skinpos, skintri) == 1;
        inskull = surface_inside(dippos, headmodel.bnd(skull).pos, headmodel.bnd(skull).tri) == 1;
        inside  = inskin & ~inskull;
      otherwise
        ft_error('unsupported surface "%s"', surface);
    end

  case {'simbio', 'duneuro'}
    % this is a model with hexahedrons or tetrahedrons
    if isfield(headmodel, 'tet')
      % the subsequent code works both for tetrahedrons or hexahedrons, but assumes the volume elements to be called "hex"
      headmodel.hex = headmodel.tet;
      headmodel = rmfield(headmodel, 'tet');
    end

    % determine the size of the relevant elements
    numhex = size(headmodel.hex,1);
    numpos = size(headmodel.pos,1);
    numdip = size(dippos,1);

    % select only the gray matter or brain tissues
    brain = find(ismember(headmodel.tissuelabel, {'gm', 'gray', 'grey', 'brain'}));

    for i=1:numel(brain)
      fprintf('flagging dipoles with tissue type "%s" as inside\n', headmodel.tissuelabel{brain(i)});
    end

    % determine all hexahedrons that are labeled as gray matter or brain tissues
    insidehex = ismember(headmodel.tissue, brain);

    % prune the mesh, only retain hexahedrons labeled as gray matter or brain tissues
    fprintf('pruning headmodel volume elements from %d to %d (%d%%)\n', numhex, sum(insidehex), round(100*sum(insidehex)/numhex));
    headmodel.hex    = headmodel.hex(insidehex,:);
    headmodel.tissue = headmodel.tissue(insidehex);
    numhex = sum(insidehex);

    % remove these, we don't need them any more
    clear brain insidehex

    % prune the mesh, i.e. only retain vertices that are part of a hexahedron
    [headmodel.pos, headmodel.hex] = remove_unused_vertices(headmodel.pos, headmodel.hex);

    % construct a sparse matrix with the mapping between all hexahedrons and vertices
    i = repmat(transpose(1:numhex), 1, size(headmodel.hex,2));
    j = headmodel.hex;
    s = ones(size(i));
    hex2pos = sparse(i(:),j(:),s(:),numhex,numpos);

    % determine the bounding box
    minpos = min(headmodel.pos,[],1);
    maxpos = max(headmodel.pos,[],1);
    insidedip = all(bsxfun(@ge, dippos, minpos),2) & all(bsxfun(@le, dippos, maxpos),2);
    fprintf('pruning dipole positions from %d to %d (%d%%)\n', numdip, sum(insidedip), round(100*sum(insidedip)/numdip));
    insidedip  = find(insidedip);
    dippos = dippos(insidedip,:);

    % find the nearest vertex for each of the dipoles in the headmodel mesh
    posindx = knnsearch(headmodel.pos, dippos);

    % The following code is only guaranteed to work with convex elements. Regular
    % hexahedra and tetrahedra are convex, and the adapted hexahedra we can use with
    % SIMBIO/Duneuro have to be convex as well.

    inside = false(1, numdip);
    % for each dipole determine whether it is inside one of the neighbouring hexahedrons
    % this will be the case for all vertices that are inside the middle, but not at the edges
    for i=1:numel(insidedip)
      hexindx = find(hex2pos(:,posindx(i)));
      for j=1:numel(hexindx)
        posx = headmodel.pos(headmodel.hex(hexindx(j),:),1);
        posy = headmodel.pos(headmodel.hex(hexindx(j),:),2);
        posz = headmodel.pos(headmodel.hex(hexindx(j),:),3);
        pos = [posx posy posz];
        if isequal(convhull(pos), convhull([pos; dippos(i,:)]))
          inside(insidedip(i)) = true;
          break % out of the for-loop
        end % if
      end % for each hexahedron
    end % for each of the dipole positions

  otherwise
    ft_error('unrecognized volume conductor model');
end

% ensure that it is a boolean column vector
inside(isnan(inside(:))) = 0;
inside = logical(inside(:));
