function sourcemodel = hartmut_add_eyes(sourcemodel, cfg, sens, headmodel, leadfieldopt)

% HARTMUT_ADD_EYES adds ocular source candidates to the source model that is used for the
% HArtMuT grid search in FT_DIPOLEFITTING. Each ocular candidate represents a mirror-symmetric
% dipole pair, and its fused leadfield L(p) + L(mirror(p))*mommap is precomputed here so that
% the grid search can treat it as an ordinary source. The existing brain and scalp leadfields
% are left empty, so they are computed on the fly during the scan.
%
% Use as
%   sourcemodel = hartmut_add_eyes(sourcemodel, cfg, sens, headmodel, leadfieldopt)
%
% When the eye position or the symmetry axis cannot be determined, the source model is
% returned unchanged, so the grid search falls back to brain and scalp candidates only.
%
% See also FT_DIPOLEFITTING, HARTMUT_EYEMODEL, FT_COMPUTE_LEADFIELD

% Copyright (C) 2026, Nils Harmening
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

% determine the candidate positions in one eye, the left-right symmetry axis, and the eye centre
coordsys = ft_getopt(headmodel, 'coordsys', ft_getopt(sens, 'coordsys', 'unknown'));
[eyepos, lraxis, centre, radius] = hartmut_eyemodel(cfg.dipfit.eye, coordsys, sourcemodel.unit);
if isempty(eyepos)
  % the eye position or symmetry axis could not be determined, keep brain and scalp only
  ft_warning('the HArtMuT ocular source positions could not be determined, the eyes are fit as single dipoles; provide cfg.dipfit.eye.pos or use an MNI-like coordinate system');
  return
end

% the partner dipole is the mirror image across the midsagittal plane, i.e. the plane where
% the left-right coordinate is zero
mommap   = ft_getopt(cfg.dipfit.constr, 'mommap', eye(3));
axisindx = find(strcmp({'x', 'y', 'z'}, lraxis));
mirror   = ones(1,3);
mirror(axisindx) = -1;

% precompute the fused leadfield for each ocular candidate
neye  = size(eyepos,1);
eyelf = cell(neye,1);
for i=1:neye
  lf1 = ft_compute_leadfield(eyepos(i,:),         sens, headmodel, leadfieldopt{:});
  lf2 = ft_compute_leadfield(eyepos(i,:).*mirror, sens, headmodel, leadfieldopt{:});
  eyelf{i} = lf1 + lf2*mommap;
end

% leave the existing brain and scalp leadfields empty, they are computed during the scan
if ~isfield(sourcemodel, 'leadfield')
  sourcemodel.leadfield = cell(size(sourcemodel.pos,1), 1);
end
if ~isfield(sourcemodel, 'compartment')
  sourcemodel.compartment = repmat({''}, size(sourcemodel.pos,1), 1);
end
sourcemodel.inside = sourcemodel.inside(:);

% remove the existing single-dipole candidates inside either eye, since the eye region is
% now represented by the symmetric pairs and would otherwise be seeded by the wrong model
d2left  = sum((sourcemodel.pos - centre).^2, 2);
d2right = sum((sourcemodel.pos - centre.*mirror).^2, 2);
ineye   = d2left<=radius^2 | d2right<=radius^2;
sourcemodel.pos(ineye,:)       = [];
sourcemodel.inside(ineye)      = [];
sourcemodel.leadfield(ineye)   = [];
sourcemodel.compartment(ineye) = [];

% merge the ocular candidates into the source model
sourcemodel.pos         = [sourcemodel.pos;         eyepos];
sourcemodel.inside      = [sourcemodel.inside;      true(neye,1)];
sourcemodel.leadfield   = [sourcemodel.leadfield;   eyelf];
sourcemodel.compartment = [sourcemodel.compartment; repmat({'eye'}, neye, 1)];
if isfield(sourcemodel, 'dim')
  sourcemodel = rmfield(sourcemodel, 'dim'); % the positions are no longer a regular grid
end
