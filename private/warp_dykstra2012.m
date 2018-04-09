function [coord_snapped] = warp_dykstra2012(cfg, elec, surf)

% WARP_DYKSTRA2012 projects the ECoG grid / strip onto a cortex hull
% using the algorithm described in Dykstra et al. (2012,
% Neuroimage) in which the distance from original positions and the
% deformation of the grid are minimized. This function relies on MATLAB's 
% optimization toolbox. To align ECoG electrodes to the pial surface, you 
% first need to compute the cortex hull with FT_PREPARE_MESH.
%
% See also FT_ELECTRODEREALIGN, FT_PREPARE_MESH, WARP_HERMES2010

% Copyright (C) 2012-2017, Andrew Dykstra, Gio Piantoni, Arjen Stolk
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

if strcmp(cfg.pairmethod, 'label') % FIXME: it were cleaner if this were housed under the create_elecpairs subfunction
  % determine if any electrodes are missing and
  % re-order electrode positions based on numbers in the label and replace
  % missing numbers with electrodes at position [NaN NaN NaN]
  digits = regexp(elec.label, '\d+', 'match');
  elec.maxdigit = 1;
  for l=1:numel(digits)
    ElecStrs{l,1} = regexprep(elec.label{l}, '\d+(?:_(?=\d))?', ''); % without electrode numbers
    labels{l,1} = digits{l}{1}; % use first found digit
    if str2num(digits{l}{1}) > elec.maxdigit
      elec.maxdigit = str2num(digits{l}{1});
    end
  end
  elec.ElecStr = cell2mat(unique(ElecStrs));
  
  pos_ordered = [];
  labels_ordered = {};
  elec.cutout = []; % index of electrodes that appear to be cut out
  dowarn = 0;
  for e = 1:elec.maxdigit
    labels_ordered{e} = [elec.ElecStr num2str(e)];
    if ~isempty(match_str(labels, num2str(e)))
      pos_ordered(e, :) = elec.elecpos(match_str(labels, num2str(e)),:);
    else
      elec.cutout(end+1) = e;
      pos_ordered(e, :) = NaN(1,3);
      dowarn = 1;
    end
  end
  if dowarn
    ft_warning('%s appears to be missing electrodes %s or have electrodes that are labeled using an unconventional numbering system', elec.ElecStr, num2str(elec.cutout));
  end
  elec.label = labels_ordered;
  elec.elecpos = pos_ordered;
end

% compute pairs of neighbors
pairs = create_elecpairs(elec, cfg.pairmethod);

if strcmp(cfg.pairmethod, 'label') % FIXME: it were cleaner if this were housed under the create_elecpairs subfunction
  % remove electrodes that are cut out from elec.elecpos and update the
  % pairs list to reflect these removals
  elec.elecpos(elec.cutout, :) = [];
  for c = length(elec.cutout):-1:1 % for each of the cutout electrodes, starting with the highest numbered one
    % find pairs that refered to elec.cutout(c) and remove them
    for n = size(pairs, 1):-1:1 % for each of the rows in pairs, starting with the highest
      if any(intersect(pairs(n, :), elec.cutout(c)))
        pairs(n, :) = [];
      end
    end
    
    % subtract 1 from all the numbers in pairs that are higher than elec.cutout(c)
    pairs(find(pairs > elec.cutout(c))) = pairs(find(pairs > elec.cutout(c)))-1;
  end
  
  % remove any pairs that reference an electrode higher than the number of
  % electrodes in elec.elecpos
  for n = size(pairs, 1):-1:1 % for each of the rows in pairs, starting with the highest
    if any(pairs(n, :) > size(elec.elecpos, 1))
      pairs(n, :) = [];
    end
  end
end


% get starting coordinates
coord0 = elec.elecpos;
coord = elec.elecpos;

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

if strcmp(cfg.pairmethod, 'label')
  % return the order of the coordinates to its original order (before they
  % were ordered sequentially)
  
  % first add back NaN's for the cutout electrodes
  tmp = coord_snapped;
  for c = 1:numel(elec.cutout)
    if elec.cutout(c) < elec.maxdigit
      tmp(elec.cutout(c)+1:end+1, :) = tmp(elec.cutout(c):end, :);
      tmp(elec.cutout(c), :) = NaN(1,3);
    end
  end
  
  for n = 1:size(coord_snapped, 1)
    coord_snapped(n, :) = tmp(str2num(cell2mat(labels(n))),:);
  end
end

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

if strcmp(method, 'label')
  
  % determine grid dimensions (1st dim: number of arrays, 2nd dim: number of elecs in an array)
  fprintf('creating electrode pairs based on electrode labels\n');
  GridDim = determine_griddim(elec);
  if GridDim(1)*GridDim(2) ~= elec.maxdigit
    warning('the product of the dimensions does not equal the maximum digit in the electrode labels, so if incorrect, use cfg.pairmethod = ''pos'' instead\n');
  elseif any(GridDim(:)==1) % if not because of strips, this could happen in case of missing electrodes
    warning('if this not a strip, there may be electrodes missing, so if incorrect, use cfg.pairmethod = ''pos'' instead\n');
  end
  
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
 
elseif strcmp(method, 'pos')
  
  % KNN_PAIRS compute pairs of neighbors of the grid
  fprintf('creating electrode pairs based on electrode positions\n');
  pairs = knn_pairs(elec.elecpos, 4);
  
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
