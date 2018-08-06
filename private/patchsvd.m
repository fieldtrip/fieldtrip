function [grid] = patchsvd(cfg, grid)

% This function does something
% cf. Limpiti et al IEEE trans biomed eng 2006;53(9);1740-54

% Copyright (c) 2006, Jan-Mathijs Schoffelen & Robert Oostenveld, F.C. Donders Centre
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

% set the defaults
if ~isfield(cfg, 'patchsvd'),    cfg.patchsvd    = 3;         end
if ~isfield(cfg, 'patchsvdnum'), cfg.patchsvdnum = 5;         end
if ~isfield(cfg, 'feedback'),    cfg.feedback    = 'textbar'; end

if isnumeric(cfg.patchsvd),
  Ndipoles = size(grid.pos,1);
else 
  Ndipoles = 1;
end
Ninside  = length(grid.inside);
Nchans   = size(grid.leadfield{grid.inside(1)}, 1);
lfall    = cell(1,Ndipoles);
coeff    = cell(1,Ndipoles);
nghbr    = cell(1,Ndipoles);

if isnumeric(cfg.patchsvd) && ~isfield(grid, 'patchindx'),
  fprintf('computing patches in 3D, not taking topology into account\n');
  progress('init', cfg.feedback, 'computing patchsvd');
  for dipindx=1:Ninside
    % renumber the loop-index variable to make it easier to print the progress bar
    i = grid.inside(dipindx);
    
    % compute the distance from this dipole to each other dipole
    dist = sqrt(sum((grid.pos-repmat(grid.pos(i,:), [Ndipoles 1])).^2, 2));
  
    % define the region of interest around this dipole
    sel  = find(dist<=cfg.patchsvd);
    sel  = intersect(sel, grid.inside);
    Nsel = length(sel);
    
    progress(dipindx/Ninside, 'computing patchsvd %d/%d, Nsel=%d\n', dipindx, Ninside, Nsel);
    % concatenate the leadfield of all dipoles that are inside the ROI into one matrix
    lfr     = cell2mat(grid.leadfield(sel(:)'));
    % svd of leadfields of dipoles inside the ROI
    [U,S,V] = svd(lfr);
    
    if cfg.patchsvdnum < 1, 
      % Limpiti et al 2006 formula 12
      s          = diag(S).^2;
      s          = cumsum(s)./sum(s);
      n(dipindx) = find(s - cfg.patchsvdnum > 0, 1);
    else
      n(dipindx) = cfg.patchsvdnum;
    end
    lfall{i} = U(:,1:n(dipindx));%*S(1:n(dipindx),1:n(dipindx)); %klopt dit?
    nghbr{i} = sel;
    coeff{i} = V(1:n(dipindx), :);
  end
  progress('close');
elseif isnumeric(cfg.patchsvd) && isfield(grid, 'patchindx'),
  fprintf('computing patches in 2D, taking surface topology into account\n');
  progress('init', cfg.feedback, 'computing patchsvd');
  for dipindx=1:Ninside
    % renumber the loop-index variable to make it easier to print the progress bar
    i = grid.inside(dipindx);
    
    % define the region of interest around this dipole
    sel  = nearest(grid.patchsize(i,:), cfg.patchsvd);
    sel  = grid.patchindx(i,1:sel);
    Nsel = length(sel);
    
    progress(dipindx/Ninside, 'computing patchsvd %d/%d, Nsel=%d\n', dipindx, Ninside, Nsel);
    % concatenate the leadfield of all dipoles that are inside the ROI into one matrix
    lfr     = cell2mat(grid.leadfield(sel(:)'));
    % svd of leadfields of dipoles inside the ROI
    [U,S,V] = svd(lfr);
    
    if cfg.patchsvdnum < 1, 
      % Limpiti et al 2006 formula 12
      s          = diag(S).^2;
      s          = cumsum(s)./sum(s);
      n(dipindx) = find(s - cfg.patchsvdnum > 0, 1);
    else
      n(dipindx) = min(cfg.patchsvdnum, size(lfr,2));
    end
    lfall{i} = U(:,1:n(dipindx));%*S(1:n(dipindx),1:n(dipindx)); %klopt dit?
    nghbr{i} = sel;
    coeff{i} = V(1:n(dipindx), :);
    sv{i}    = s;
  end
  progress('close');
elseif strcmp(cfg.patchsvd, 'all'),
  lfr     = cell2mat(grid.leadfield(grid.inside(:)'));
  [U,S,V] = svd(lfr);
  
  if cfg.patchsvdnum < 1, 
    % Limpiti et al 2006 formula 12
    s = diag(S).^2;
    s = cumsum(s)./sum(s);
    n = find(s - cfg.patchsvdnum > 0, 1);
  else
    n = cfg.patchsvdnum;
  end
  lfall{1} = U(:,1:n);%*S(1:n(dipindx),1:n(dipindx)); %klopt dit?
  nghbr{1} = grid.inside;
  coeff{1} = V(1:n, :);
 
  %---change output
  grid.pos    = mean(grid.pos(grid.inside,:),1);
  grid.inside = 1; 
  grid.dim    = [1 1 1];
  grid.xgrid  = grid.pos(1);
  grid.ygrid  = grid.pos(2);
  grid.zgrid  = grid.pos(3); 
else
  %do nothing
end

if ~all(n==n(1)),
  nmax = max(n);
  for dipindx = 1:Ninside
    i = grid.inside(dipindx);
    if n(dipindx)<nmax,
      lfall{i} = [lfall{i} zeros(Nchans, nmax-n(dipindx))];
    end
  end
end

% update the leadfields
grid.leadfield = lfall;
%grid.expcoeff  = coeff;
grid.neighbours= nghbr;
