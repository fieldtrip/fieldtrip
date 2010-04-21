function bnd = prepare_mesh_headshape(cfg)

% PREPARE_MESH_HEADSHAPE
%
% See also PREPARE_MESH_MANUAL, PREPARE_MESH_SEGMENTATION

% Copyrights (C) 2009, Robert Oostenveld
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

% get the surface describing the head shape
if isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
  % use the headshape surface specified in the configuration
  headshape = cfg.headshape;
elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
  % use the headshape points specified in the configuration
  headshape.pnt = cfg.headshape;
elseif ischar(cfg.headshape)
  % read the headshape from file
  headshape = read_headshape(cfg.headshape);
else
  error('cfg.headshape is not specified correctly')
end

% usually a headshape only describes a single surface boundaries, but there are cases
% that multiple surfaces are included, e.g. skin, skull, brain
nbnd = numel(headshape);

if ~isfield(headshape, 'tri')
  % generate a closed triangulation from the surface points
  for i=1:nbnd
    headshape(i).pnt = unique(headshape(i).pnt, 'rows');
    headshape(i).tri = projecttri(headshape(i).pnt);
  end
end

if ~isempty(cfg.numvertices) && ~strcmp(cfg.numvertices, 'same')
  for i=1:nbnd
    [tri1, pnt1] = reducepatch(headshape(i).tri, headshape(i).pnt, 3*cfg.numvertices);
    % remove double vertices
    pnt1 = unique(pnt1, 'rows');
    % reconstruct the triangulation
    tri1 = projecttri(pnt1);
    % replace the probably unevenly distributed triangulation with a regular one
    % and retriangulate it to the desired accuracy
    [pnt2, tri2] = msphere(cfg.numvertices); % this is a regular triangulation
    [headshape(i).pnt, headshape(i).tri] = retriangulate(pnt1, tri1, pnt2, tri2, 2);
  end
end

% the output should only describe one or multiple boundaries and should not
% include any other fields
bnd = rmfield(headshape, setdiff(fieldnames(headshape), {'pnt', 'tri'}));
