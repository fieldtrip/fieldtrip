function [coord_snapped] = warp_dykstra2012(cfg, elec, surf)

% WARP_DYKSTRA2012 projects the ECoG grid / strip onto a cortex hull
% while minimizing the distance from original positions and the
% deformation of the grid. To align ECoG electrodes to the pial surface,
% you first need to compute the cortex hull with FT_PREPARE_MESH.
%
% WARP_DYKSTRA2012 uses the algorithm described in Dykstra et al. (2012,
% Neuroimage) in which electrodes are projected onto pial surface while
% minimizing the displacement of the electrodes from original location
% and maintaining the grid shape. It relies on the optimization toolbox.
%
% See also FT_ELECTRODEREALIGN, FT_PREPARE_MESH

% Copyright (C) 2012-2017, Arjen Stolk, Gio Piantoni, Andrew Dykstra
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

disp('using warp algorithm described in Dykstra et al. 2012 Neuroimage PMID: 22155045')

% set the defaults
cfg.feedback      = ft_getopt(cfg, 'feedback', 'no');

% undocumented local options
cfg.pairmethod    = ft_getopt(cfg, 'pairmethod', 'pos'); % eletrode pairing based on electrode 'pos' or 'label' (for computing deformation energy)
cfg.deformweight  = ft_getopt(cfg, 'deformweight',  1); % weight of deformation relative to shift energy cost

% get starting coordinates
coord0 = elec.elecpos;
coord = elec.elecpos;

% compute pairs of neighbors
pairs = create_elecpairs(elec, cfg.pairmethod);

% anonymous function handles
efun = @(coord_snapped) energy_electrodesnap(coord_snapped, coord, pairs, cfg.deformweight);
cfun = @(coord_snapped) dist_to_surface(coord_snapped, surf);

% options
%   'UseParallel', 'always',...
options = optimset('Algorithm','active-set',...
  'MaxIter', 50,...
  'MaxFunEvals', Inf,...
  'GradObj', 'off',...
  'TypicalX', coord(:),...
  'DiffMaxChange', 2,...
  'DiffMinChange', 0.3,...
  'TolFun', 0.3,...
  'TolCon', 0.01 * size(coord0, 1),...
  'TolX', 0.5,...
  'Diagnostics', 'off',...
  'RelLineSrchBnd',1);

if strcmp(cfg.feedback, 'yes')
  options = optimset(options, 'Display', 'iter');
else
  options = optimset(options, 'Display', 'final');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy Minimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run minimization: efun (shift + deform energy) is minimized; cfun (surface distance) is a nonlinear constraint
coord_snapped = fmincon(efun, coord0, [], [], [], [], [], [], cfun, options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function energy = energy_electrodesnap(coord, coord_orig, pairs, deformweight) % (minimized) energy function
% ENERGY_ELECTRODESNAP compute energy to be minimized, based on deformation
% and distance of the electrodes from original distance

% energy needed to move electrodes to the surface
energy_eshift = sum((coord - coord_orig).^2, 2);

% energy needed to deform grid shape
dist = sqrt(sum((coord(pairs(:, 1), :) - coord(pairs(:, 2), :)).^2, 2));
dist_orig = sqrt(sum((coord_orig(pairs(:, 1), :) - coord_orig(pairs(:, 2), :)).^2, 2));
energy_deform = (dist - dist_orig).^2;

% (weighted) sum of the above
energy = mean(energy_eshift) + (deformweight * mean(energy_deform.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, dist] = dist_to_surface(coord, surf) % (nonlinear) constraint function
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pairs = create_elecpairs(elec, method)

if strcmp(method, 'label'); 
  
  % determine grid dimensions (1st dim: number of arrays, 2nd dim: number of elecs in an array)
  fprintf('creating electrode pairs based on electrode labels\n');
  GridDim = determine_griddim(elec);
  
  % create pairs based on dimensions
  diagonal = 1;
  pairs = [];
  for e = 1:GridDim(1)*GridDim(2)   
    if isequal(mod(e,GridDim(2)), 1) % begin of each elec array
      pairs(end+1,:) = [e e+1]; % following elec
      pairs(end+1,:) = [e e-GridDim(2)]; % adjacent preceding elec
      pairs(end+1,:) = [e e+GridDim(2)]; % adjacent following elec
      if diagonal
        pairs(end+1,:) = [e e-GridDim(2)+1]; % adjacent preceding elec
        pairs(end+1,:) = [e e+GridDim(2)+1]; % adjacent following elec
      end
    elseif isequal(mod(e,GridDim(2)), 0) % end of each elec array
      pairs(end+1,:) = [e e-1]; % preceding elec
      pairs(end+1,:) = [e e-GridDim(2)]; % adjacent preceding elec
      pairs(end+1,:) = [e e+GridDim(2)]; % adjacent following elec
      if diagonal
        pairs(end+1,:) = [e e-GridDim(2)-1]; % adjacent preceding elec
        pairs(end+1,:) = [e e+GridDim(2)-1]; % adjacent following elec
      end
    else
      pairs(end+1,:) = [e e-1]; % preceding elec
      pairs(end+1,:) = [e e+1]; % following elec
      pairs(end+1,:) = [e e-GridDim(2)]; % adjacent preceding elec
      pairs(end+1,:) = [e e+GridDim(2)]; % adjacent following elec
      if diagonal
        pairs(end+1,:) = [e e-GridDim(2)-1]; % adjacent preceding elec
        pairs(end+1,:) = [e e+GridDim(2)-1]; % adjacent following elec
        pairs(end+1,:) = [e e-GridDim(2)+1]; % adjacent preceding elec
        pairs(end+1,:) = [e e+GridDim(2)+1]; % adjacent following elec
      end
    end
  end
  pairs( pairs(:,2)<1 | pairs(:,2)>GridDim(1)*GridDim(2) ,:) = []; % out of bounds
  pairs = unique(sort(pairs,2),'rows'); % unique pairs
 
elseif strcmp(method, 'pos');
  
  % KNN_PAIRS compute pairs of neighbors of the grid
  fprintf('creating electrode pairs based on electrode positions\n');
  pairs = knn_pairs(elec.elecpos, 4);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GridDim = determine_griddim(elec)
% assumes intact grids/strips and elec count starting at 1
% A. Stolk, 2017

% extract numbers from elec labels
digits = regexp(elec.label, '\d+', 'match');
for l=1:numel(digits)
  labels{l,1} = digits{l}{1}; % use first found digit
end

% determine grid dimensions (1st dim: number of arrays, 2nd dim: number of elecs in an array)
if isequal(numel(labels), 256)
  GridDim(1) = 16; GridDim(2) = 16;
elseif isequal(numel(labels), 64)
  GridDim(1) = 8; GridDim(2) = 8;
elseif isequal(numel(labels), 48)
  e6 = elec.elecpos(match_str(labels, num2str(6)),:);
  e7 = elec.elecpos(match_str(labels, num2str(7)),:);
  e8 = elec.elecpos(match_str(labels, num2str(8)),:);
  e9 = elec.elecpos(match_str(labels, num2str(9)),:);
  d6to7 = sqrt(sum((e6-e7).^2)); % distance of elec 6 to 7
  d8to9 = sqrt(sum((e8-e9).^2)); % distance of elec 8 to 9
  if d8to9 >= d6to7  % break between e8 and e9
    GridDim(1) = 6; GridDim(2) = 8;
  elseif d6to7 > d8to9 % break between e6 and e7
    GridDim(1) = 8; GridDim(2) = 6;
  end
elseif isequal(numel(labels), 32)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e8 = elec.elecpos(match_str(labels, num2str(8)),:);
  e9 = elec.elecpos(match_str(labels, num2str(9)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d8to9 = sqrt(sum((e8-e9).^2)); % distance of elec 8 to 9
  if d8to9 >= d4to5  % break between e8 and e9
    GridDim(1) = 4; GridDim(2) = 8;
  elseif d4to5 > d8to9 % break between e4 and e5
    GridDim(1) = 8; GridDim(2) = 4;
  end
elseif isequal(numel(labels), 24)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e6 = elec.elecpos(match_str(labels, num2str(6)),:);
  e7 = elec.elecpos(match_str(labels, num2str(7)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d6to7 = sqrt(sum((e6-e7).^2)); % distance of elec 6 to 7
  if d6to7 >= d4to5 % break between e6 and e7
    GridDim(1) = 4; GridDim(2) = 6;
  elseif d4to5 > d6to7 % break between e4 and e5
    GridDim(1) = 6; GridDim(2) = 4;
  end
elseif isequal(numel(labels), 20)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e6 = elec.elecpos(match_str(labels, num2str(6)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d5to6 = sqrt(sum((e5-e6).^2)); % distance of elec 5 to 6
  if d5to6 >= d4to5 % break between e5 and e6
    GridDim(1) = 4; GridDim(2) = 5;
  elseif d4to5 > d5to6 % break between e4 and e5
    GridDim(1) = 5; GridDim(2) = 4;
  end
elseif isequal(numel(labels), 16)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e8 = elec.elecpos(match_str(labels, num2str(8)),:);
  e9 = elec.elecpos(match_str(labels, num2str(9)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d8to9 = sqrt(sum((e8-e9).^2)); % distance of elec 8 to 9
  d4to8 = sqrt(sum((e4-e8).^2)); % distance of elec 4 to 8
  if d8to9 > 2*d4to5 % break between e8 and e9
    GridDim(1) = 2; GridDim(2) = 8;
  elseif d4to5 > 2*d4to8 % break between e4 and e5
    GridDim(1) = 4; GridDim(2) = 4;
  else
    GridDim(1) = 1; GridDim(2) = 16;
  end
elseif isequal(numel(labels), 12)
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  e6 = elec.elecpos(match_str(labels, num2str(6)),:);
  e7 = elec.elecpos(match_str(labels, num2str(7)),:);
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  d5to6 = sqrt(sum((e5-e6).^2)); % distance of elec 5 to 6
  d6to7 = sqrt(sum((e6-e7).^2)); % distance of elec 6 to 7
  if d4to5 > 2*d5to6 % break between e4 and e5
    GridDim(1) = 3; GridDim(2) = 4; % 4x3 unsuppported
  elseif d6to7 > 2*d5to6 % break between e6 and e7
    GridDim(1) = 2; GridDim(2) = 6;
  else
    GridDim(1) = 1; GridDim(2) = 12;
  end
elseif isequal(numel(labels), 8)
  e3 = elec.elecpos(match_str(labels, num2str(3)),:);
  e4 = elec.elecpos(match_str(labels, num2str(4)),:);
  e5 = elec.elecpos(match_str(labels, num2str(5)),:);
  d3to4 = sqrt(sum((e3-e4).^2)); % distance of elec 3 to 4
  d4to5 = sqrt(sum((e4-e5).^2)); % distance of elec 4 to 5
  if d4to5 > 2*d3to4 % break between e4 and e5
    GridDim(1) = 2; GridDim(2) = 4;
  else
    GridDim(1) = 1; GridDim(2) = 8;
  end
else
  GridDim(1) = 1; GridDim(2) = numel(labels);
end

% provide feedback of what grid dimensions were found
if any(GridDim==1) % if not because of strips, this could happen in case of missing electrodes
  warning('assuming %d x %d grid dimensions: if incorrect, use cfg.pairmethod = ''pos'' instead\n', GridDim(1), GridDim(2));
else
  fprintf('assuming %d x %d grid dimensions\n', GridDim(1), GridDim(2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pairs = knn_pairs(coord, k)
% KNN_PAIRS compute pairs of neighbors of the grid

knn_ind = knn_search(coord, coord, k);
pairs = cat(3, knn_ind, repmat([1:size(coord,1)]',1,k));
pairs = permute(pairs,[3 1 2]);
pairs = sort(reshape(pairs,2,[]),1)';
pairs = unique(pairs,'rows');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
