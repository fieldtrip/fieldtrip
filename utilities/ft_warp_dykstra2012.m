function [coord_snapped] = ft_warp_dykstra2012(coord, surf, feedback)

% FT_WARP_DYKSTRA2012 projects the ECoG grid / strip onto a cortex hull
% while minimizing the distance from original positions and the
% deformation of the grid. To align ECoG electrodes to the pial surface,
% you first need to compute the cortex hull with FT_PREPARE_MESH.
% FT_WARP_DYKSTRA2012 uses algorithm described in Dykstra et al. (2012,
% Neuroimage) in which electrodes are projected onto pial surface while
% minimizing the displacement of the electrodes from original location
% and maintaining the grid shape. It relies on the optimization toolbox.
%
% See also FT_ELECTRODEREALIGN, FT_PREPARE_MESH

% Copyright (C) 2012-2016, Gio Piantoni, Andrew Dykstra
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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

% determine whether the MATLAB Optimization toolbox is available and can be used
ft_hastoolbox('optim', 1);

disp('Please cite: Dykstra et al. 2012 Neuroimage PMID: 22155045')

% get starting coordinates
coord0 = coord;

% compute pairs of neighbors
pairs = knn_pairs(coord, 4);

% anonymous function handles
efun = @(coord_snapped) energy_electrodesnap(coord_snapped, coord, pairs);
cfun = @(coord_snapped) dist_to_surface(coord_snapped, surf);

% options
options = optimset('Algorithm','active-set',...
  'MaxIter', 50,...
  'MaxFunEvals', Inf,...
  'UseParallel', 'always',...
  'GradObj', 'off',...
  'TypicalX', coord(:),...
  'DiffMaxChange', 2,...
  'DiffMinChange', 0.3,...
  'TolFun', 0.3,...
  'TolCon', 0.01 * size(coord0, 1),...
  'TolX', 0.5,...
  'Diagnostics', 'off',...
  'RelLineSrchBnd',1);

if strcmp(feedback, 'yes')
  options = optimset(options, 'Display', 'iter');
else
  options = optimset(options, 'Display', 'final');
end

% run minimization
coord_snapped = fmincon(efun, coord0, [], [], [], [], [], [], cfun, options);

end

function [energy, denergy] = energy_electrodesnap(coord, coord_orig, pairs)
% ENERGY_ELECTRODESNAP compute energy to be minimized, based on deformation
% and distance of the electrodes from original distance

energy_eshift = sum((coord - coord_orig).^2, 2);
energy_deform = deformation_energy(coord, coord_orig, pairs);
energy = mean(energy_eshift) + mean(energy_deform.^2);

denergy=[];

end

function energy = deformation_energy(coord, coord_orig, pairs)
% DEFORMATION_ENERGY measure energy due to grid deformation

dist = sqrt(sum((coord(pairs(:, 1), :) - coord(pairs(:, 2), :)) .^ 2, 2));
dist_orig = sqrt(sum((coord_orig(pairs(:, 1), :) - coord_orig(pairs(:, 2), :)) .^ 2, 2));

energy = (dist - dist_orig) .^2;
end

function [c, dist] = dist_to_surface(coord, surf)
% DIST_TO_SURFACE Compute distance to surface, this is the fastest way to run
% it, although running the loops in other directions might be more intuitive.

c = [];

dist = zeros(size(coord, 1), 1);
for i0 = 1:size(coord, 1)
  dist_one_elec = zeros(size(surf.pos, 1), 1);
  for i1 = 1:size(surf.pos, 2)
    dist_one_elec = dist_one_elec + (surf.pos(:, i1) - coord(i0, i1)) .^ 2;
  end
  dist(i0) = min(dist_one_elec);
end
dist = sqrt(dist);

end

function pairs = knn_pairs(coord, k)
% KNN_PAIRS compute pairs of neighbors of the grid

knn_ind = knn_search(coord, coord, k);
pairs = cat(3, knn_ind, repmat([1:size(coord,1)]',1,k));
pairs = permute(pairs,[3 1 2]);
pairs = sort(reshape(pairs,2,[]),1)';
pairs = unique(pairs,'rows');

end

function idx = knn_search(Q, R, K)
%KNN_SEARCH perform search using k-Nearest Neighbors algorithm

[N, M] = size(Q);
L = size(R, 1);
idx = zeros(N, K);
D = idx;

for k = 1:N
  d = zeros(L, 1);
  for t = 1:M
    d = d + (R(:, t) - Q(k, t)) .^ 2;
  end
  d(k) = inf;
  [s, t] = sort(d);
  idx(k, :) = t(1:K);
  D(k, :)= s(1:K);
end

end
