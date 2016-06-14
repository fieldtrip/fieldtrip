function [coord_snapped] = ft_warp_dykstra2012(coord, surf, feedback)

% FT_WARP_OPTIM determine intermediate positions using warping (deformation)
% the input cloud of points is warped to match the target.
% The strategy is to start with simpelest linear warp, followed by a more
% elaborate linear warp, which then is followed by the nonlinear warps up
% to the desired order.
%
% [result, M] = ft_warp_pnt(input, target, method)
%     input          contains the Nx3 measured 3D positions
%     target         contains the Nx3 template 3D positions
%     method         should be any of 
%                     'rigidbody'
%                     'globalrescale'
%                     'traditional' (default)
%                     'nonlin1'
%                     'nonlin2'
%                     'nonlin3'
%                     'nonlin4'
%                     'nonlin5'
%
% The default warping method is a traditional linear warp with individual
% rescaling in each dimension. Optionally you can select a nonlinear warp
% of the 1st (affine) up to the 5th order.
%
% When available, this function will use the MATLAB optimization toolbox.
%
% See also FT_WARP_APPLY, FT_WARP_ERRROR

% Copyright (C) 2016, Gio Piantoni
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

energy_eshift = sum((coord - coord_orig).^2, 2);

energy_deform = deformation_energy(coord, coord_orig, pairs);

energy = mean(energy_eshift) + mean(energy_deform.^2);

denergy=[];

end

function energy = deformation_energy(coord, coord_orig, pairs)

dist = sqrt(sum((coord(pairs(:, 1), :) - coord(pairs(:, 2), :)) .^ 2, 2));
dist_orig = sqrt(sum((coord_orig(pairs(:, 1), :) - coord_orig(pairs(:, 2), :)) .^ 2, 2));

energy = (dist - dist_orig) .^2;
end

function [c, dist] = dist_to_surface(coord, surf)
% Compute distance to surface, this is the fastest way to run it, although
% running the loops in other directions might be more intuitive.

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

knn_ind = knnsearch(coord, coord, k);
pairs = cat(3, knn_ind, repmat([1:size(coord,1)]',1,k));
pairs = permute(pairs,[3 1 2]);
pairs = sort(reshape(pairs,2,[]),1)';
pairs = unique(pairs,'rows');

end

function idx = knnsearch(Q, R, K)

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

