function [stattfce, cfg] = tfcestat(cfg, statobs)

% TFCESTAT computes threshold-free cluster statistic multidimensional channel-freq-time or
% volumetric source data
%
% See also CLUSTERSTAT, FINDCLUSTER

% Copyright (C) 2021, Jan-Mathijs Schoffelen
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
cfg.feedback     = ft_getopt(cfg, 'feedback',     'text');
cfg.spmversion   = ft_getopt(cfg, 'spmversion',   'spm12');
cfg.dim          = ft_getopt(cfg, 'dim',          []);
cfg.inside       = ft_getopt(cfg, 'inside',       []);
cfg.tail         = ft_getopt(cfg, 'tail',         0);         % -1, 0, 1

cfg.tfce_h0      = ft_getopt(cfg, 'tfce_h0',      0);
cfg.tfce_H       = ft_getopt(cfg, 'tfce_H',       2);
cfg.tfce_E       = ft_getopt(cfg, 'tfce_E',       0.5);
cfg.tfce_nsteps  = ft_getopt(cfg, 'nsteps',       100);
cfg.height       = ft_getopt(cfg, 'height',       []);
% these defaults are already set in the caller function, 
% but may be necessary if a user calls this function directly
cfg.connectivity = ft_getopt(cfg, 'connectivity', false);

% ensure that the preferred SPM version is on the path
ft_hastoolbox(cfg.spmversion, 1);

if isempty(cfg.dim)
  ft_error('cfg.dim should be defined and not empty');
end

if isempty(cfg.inside)
  cfg.inside = true(cfg.dim);
end % cfg.inside is set in ft_sourcestatistics, but is also needed for timelock and freq

if isfield(cfg, 'origdim')
  cfg.dim = cfg.origdim;
end % this snippet is to support correct clustering of N-dimensional data, not fully tested yet

% get connectivity matrix for the spatially neighbouring elements
connmat = full(ft_getopt(cfg, 'connectivity', false));

needpos = cfg.tail==0 || cfg.tail== 1;
needneg = cfg.tail==0 || cfg.tail==-1;

% remove the offset, which by default is 0
statobs = statobs - cfg.tfce_h0; 

% compute the stepsize
if isempty(cfg.height)
  % explicitly compute the height here, to be added to the output cfg, so
  % that in a subsequent call to this function the same settings can be
  % used (i.e. for the randomizations)
  if needneg && needpos
    cfg.height = max(abs(statobs(:)));
  elseif needneg
    cfg.height = max(-statobs(:));
  elseif needpos
    cfg.height = max(statobs(:));
  end
end
stepsize = cfg.height./cfg.tfce_nsteps;

% first do the clustering on the observed data
spacereshapeable = (isscalar(connmat) && ~isfinite(connmat));

if needpos
  if spacereshapeable
    % this pertains to data for which the spatial dimension can be reshaped
    % into 3D, i.e. when it is described on an ordered set of positions on
    % a 3D-grid. It deals with the inside dipole positions, and creates a
    % fake extra spatial dimension, so that findcluster can deal with it
    tmp = zeros([1 cfg.dim]); 
    tmp(cfg.inside) = statobs;
  else
    tmp = reshape(statobs, [cfg.dim 1]);
  end
  
  statobspos = zeros(size(tmp));
  for j = 1:cfg.tfce_nsteps
    thr = (j-1)*stepsize;
    tmp(tmp<=thr) = 0;
    [clus, nclus] = findcluster(tmp, connmat);
    extent = zeros(size(clus));
    extent = getextent(clus, nclus, extent); 
    statobspos = statobspos + stepsize .* (extent.^cfg.tfce_E) .* (thr.^cfg.tfce_H);
  end
  
  if spacereshapeable
    statobspos = statobspos(cfg.inside);
  else
    statobspos = statobspos(:);
  end
  
end % if needpos

if needneg
  if spacereshapeable
    tmp = zeros([1 cfg.dim]); 
    tmp(cfg.inside) = statobs;
  else
    tmp = reshape(statobs, [cfg.dim 1]);
  end
  tmp = -tmp;
  
  statobsneg = zeros(size(tmp));
  for j = 1:cfg.tfce_nsteps
    thr = (j-1)*stepsize;
    tmp(tmp<=thr) = 0;
    [clus, nclus] = findcluster(tmp, connmat);
    extent = zeros(size(clus));
    extent = getextent(clus, nclus, extent);
    statobsneg = statobsneg + stepsize .* (extent.^cfg.tfce_E) .* (thr.^cfg.tfce_H);
  end
  
  if spacereshapeable
    statobsneg = -statobsneg(cfg.inside);
  else
    statobsneg = -statobsneg(:);
  end
  
end % if needneg
stattfce = statobsneg + statobspos;


% faster integration a la Bruno Giordano, borrowed from LIMO_eeg
function extent_map = getextent(clustered_map, num, extent_map)

clustered_map = clustered_map(:);
nv = histc(clustered_map,0:num);
[dum,idxall] = sort(clustered_map,'ascend');
idxall(1:nv(1)) = [];
nv(1) = [];
ends = cumsum(nv);
inis = ends-nv+1;
for i = 1:num
  idx = idxall(inis(i):ends(i));
  extent_map(idx) = nv(i);
end

