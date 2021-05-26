function [sourceout] = ft_source2full(sourcein)

% FT_SOURCE2FULL recreates the grid locations outside the brain in the source
% reconstruction, so that the source volume again describes the full grid.
% This undoes the memory savings that can be achieved using FT_SOURCE2SPARSE
% and makes it possible again to plot the source volume and save it to an
% external file.
%
% Use as
%   [source] = ft_source2full(source)
%
% See also FT_SOURCE2SPARSE

% Copyright (C) 2004-2021, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% this function does not take a cfg as input
cfg = [];

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance sourcein

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

%%

sourcein = ft_checkdata(sourcein, 'insidestyle', 'logical', 'hasdim', true);

% this assumes that the voxel data are ordered as if in a regularly spaced 3D grid,
% but with only the inside voxels present
ft_warning('assuming the positions to be ordered in a regularly spaced 3D grid');

sparsepos = sourcein.pos;
nsparse = size(sparsepos,1);

% this is a first approximation of a point in the center of the head
centerpos = mean(sparsepos,1);
dist = sparsepos - repmat(centerpos, nsparse, 1);
dist = sqrt(sum(dist.^2, 2));
[mindist, minindx] = min(dist);

% determine a single grid point in the center of the head
% it should have neighbours in all directions
centerpos = sparsepos(minindx,:);
dist = sparsepos - repmat(centerpos, nsparse, 1);
dist = sqrt(sum(dist.^2, 2));
dist = round(dist, 3, 'significant');
dist = sort(unique(dist));

dist0 = dist(1); % this is zero
dist1 = dist(2); % this is the distance over one edge
dist2 = dist(3); % this is the distance over the diagonal of a square
dist3 = dist(4); % this is the distance over the diagonal of a cube

% this assumes that the spacing between grid points is isotropic
assert(round(dist2/dist1 - sqrt(2),2)==0);
assert(round(dist3/dist1 - sqrt(3),2)==0);

dist = sparsepos - repmat(centerpos, nsparse, 1);
dist = sqrt(sum(dist.^2, 2));
dist = round(dist, 3, 'significant');

% find the nearest neighbours in each direction
nbpos = find(dist==dist1); % over one edge, there should be  6 of these
nb2 = find(dist==dist2); % over the diagonal of a square, there should be 12 of these
nb3 = find(dist==dist3); % over the diagonal of a cube, there should be  8 of these

assert(length(nbpos)==6);
assert(length(nb2)==12);
assert(length(nb3)==8);

% take the 6 nearest neighbours
nbpos = sparsepos(nbpos,:);
nbpos(:,1) = nbpos(:,1) - centerpos(1);
nbpos(:,2) = nbpos(:,2) - centerpos(2);
nbpos(:,3) = nbpos(:,3) - centerpos(3);

[mx, ix] = max(nbpos(:,1));
ux = nbpos(ix,:); ux = ux/norm(ux);
[my, iy] = max(nbpos(:,2));
uy = nbpos(iy,:); uy = uy/norm(uy);
[mz, iz] = max(nbpos(:,3));
uz = nbpos(iz,:); uz = uz/norm(uz);

u = [ux; uy; uz];

if det(u)<0
  % ensure that x, y, and z form a right-handed coordinate system
  u(3,:) = -u(3,:);
end

shiftedpos = sparsepos;
shiftedpos(:,1) = shiftedpos(:,1) - centerpos(1);
shiftedpos(:,2) = shiftedpos(:,2) - centerpos(2);
shiftedpos(:,3) = shiftedpos(:,3) - centerpos(3);

ijk = shiftedpos/u;
ijk(:,1) = ijk(:,1) - min(ijk(:,1)) + 1;
ijk(:,2) = ijk(:,2) - min(ijk(:,2)) + 1;
ijk(:,3) = ijk(:,3) - min(ijk(:,3)) + 1;
ijk = round(ijk);

dim = [max(ijk(:,1)) max(ijk(:,2)) max(ijk(:,3))];
ind = sub2ind(dim, ijk(:,1), ijk(:,2), ijk(:,3));

nfull = prod(dim);
xgrid = 1:dim(1);
ygrid = 1:dim(2);
zgrid = 1:dim(3);
[X, Y, Z] = ndgrid(xgrid,ygrid,zgrid);

gridpos = [X(:) Y(:) Z(:)];

% sparsepos' = transform * gridpos', hence transform = sparsepos' / gridpos'
% where the 4th column (or row) with ones needs to be added
transform = [sparsepos ones(nsparse,1)]' / [gridpos(ind,:) ones(nsparse,1)]';

% construct the positions of the full grid
fullpos = transform * [gridpos ones(nfull,1)]';
fullpos = fullpos(1:3,:)';

% check that the reconstructed inside positions match the original inside positions
reconstructed = fullpos(ind,:);
assert(max(abs(sparsepos(:)-reconstructed(:)))<1e-6);

%% construct the output data

% construct a boolean volume of the inside and outside grid points
inside = false(dim);
inside(ind) = true;
outside = ~inside;

% keep the descriptive fields that remain valid
sourceout = keepfields(sourcein, {'time', 'freq'});

sourceout.dim = dim;
sourceout.transform = transform;
sourceout.pos = fullpos;      % this is represented as a Nx3 matrix
sourceout.inside = inside(:); % this is represented as a Nx1 vector

% convert to vectors with indices
inside  = find(inside(:));
outside = find(outside(:));

fprintf('total number of dipoles        : %d\n', length(inside)+length(outside));
fprintf('number of dipoles inside  brain: %d\n', length(inside));
fprintf('number of dipoles outside brain: %d\n', length(outside));

% loop over parameters
fn = fieldnames(sourcein);
fn = setdiff(fn, {'dim', 'transform', 'pos', 'inside', 'cfg'});

for i=1:numel(fn)
  dimord = getdimord(sourcein, fn{i});
  
  clear tmp
  if startsWith(dimord, 'pos_pos') % this should go first
    tmp(inside,inside,:,:,:,:,:) = sourcein.(fn{i});
    tmp(outside,outside,:,:,:,:,:) = nan;
    sourceout.(fn{i}) = tmp;
  elseif startsWith(dimord, 'pos')
    tmp(inside,:,:,:,:,:,:) = sourcein.(fn{i});
    tmp(outside,:,:,:,:,:,:) = nan;
    sourceout.(fn{i}) = tmp;
  elseif startsWith(dimord, '{pos}')
    tmp = cell(nfull,1);
    tmp(inside) = sourcein.(fn{i});
    tmp(outside) = {[]}; % empty cell
    sourceout.(fn{i}) = tmp;
  elseif startsWith(dimord, '{pos_pos}')
    tmp = cell(nfull,nfull);
    tmp(inside,inside) = sourcein.(fn{i});
    tmp(outside,outside) = {[]}; % empty cell
    sourceout.(fn{i}) = tmp;
  else
    % skipping parameters with unknown dimord
  end
  
  if exist('tmp', 'var')
    ft_info('making "%s" with dimord "%s" full', fn{i}, dimord);
  else
    ft_info('skipping "%s" with dimord "%s"', fn{i}, dimord);
  end
  
end % for each field

%% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   sourcein
ft_postamble provenance sourceout
ft_postamble history    sourceout
