function [stat, cfg] = tfcestat(cfg, statrnd, statobs)

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

% these defaults are already set in the caller function, 
% but may be necessary if a user calls this function directly
cfg.connectivity     = ft_getopt(cfg, 'connectivity',     false);

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

Nrand      = size(statrnd,2);
prb_pos    = zeros(size(statobs));
prb_neg    = zeros(size(statobs));

% remove the offset, which by default is 0
statrnd = statrnd - cfg.tfce_h0;
statobs = statobs - cfg.tfce_h0; 

% compute the stepsize
if needneg && needpos
  height = max(max(abs(statobs(:))), max(abs(statrnd(:))));
elseif needneg
  height = max(max(-statobs(:)), max(-statrnd(:)));
elseif needpos
  height = max(max(statobs(:)), max(statrnd(:)));
end
stepsize = height./cfg.tfce_nsteps;

% first do the clustering on the observed data
spacereshapeable = (numel(connmat)==1 && ~isfinite(connmat));

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
    [clus, nclus] = findcluster(tmp, connmat, 0);
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
    [clus, nclus] = findcluster(tmp, connmat, 0);
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

% add the offset back to the observed data
statobs = statobs + cfg.tfce_h0;

posdistribution = zeros(1,Nrand); % this holds the statistic of the largest positive tfce value in each randomization
negdistribution = zeros(1,Nrand); % this holds the statistic of the largest negative tfce value in each randomization

% do the clustering on the randomized data
ft_progress('init', cfg.feedback, 'computing tfce for the test statistic computed from the randomized design');
for i = 1:Nrand
  ft_progress(i/Nrand, 'computing tfce for randomization %d from %d\n', i, Nrand);
  if needpos
    if spacereshapeable
      tmp = zeros([1 cfg.dim]);
      tmp(cfg.inside) = statrnd(:,i);
    else
      tmp = reshape(statrnd(:,i), [cfg.dim 1]);
    end
    
    statrndpos = zeros(size(tmp));
    for j = 1:cfg.tfce_nsteps
      thr = (j-1)*stepsize;
      tmp(tmp<=thr) = 0;
      [clus, nclus] = findcluster(tmp, connmat, 0);
      extent = zeros(size(clus));
      extent = getextent(clus, nclus, extent);
      statrndpos = statrndpos + stepsize .* (extent.^cfg.tfce_E) .* (thr.^cfg.tfce_H);
    end
    
    if spacereshapeable
      statrndpos = statrndpos(cfg.inside);
    else
      statrndpos = statrndpos(:);
    end  
    posdistribution(i) = max(statrndpos(:));
    
    prb_pos = prb_pos + (statobspos<posdistribution(i));
  end % needpos
  
  if needneg
    if spacereshapeable
      tmp = zeros([1 cfg.dim]);
      tmp(cfg.inside) = statrnd(:,i);
    else
      tmp = reshape(statrnd(:,i), [cfg.dim 1]);
    end
    tmp = -tmp;
    
    statrndneg = zeros(size(tmp));
    for j = 1:cfg.tfce_nsteps
      thr = (j-1)*stepsize;
      tmp(tmp<=thr) = 0;
      [clus, nclus] = findcluster(tmp, connmat, 0);
      extent = zeros(size(clus));
      extent = getextent(clus, nclus, extent);
      statrndneg = statrndneg + stepsize .* (extent.^cfg.tfce_E) .* (thr.^cfg.tfce_H);
    end
    
    if spacereshapeable
      statrndneg = -statrndneg(cfg.inside);
    else
      statrndneg = -statrndneg(:);
    end
    negdistribution(i) = min(statrndneg(:));
  
    prb_neg = prb_neg + (statobsneg>negdistribution(i));
  end % needneg
end % for 1:Nrand
ft_progress('close');

if isequal(cfg.numrandomization, 'all')
  N = Nrand;
else
  % the minimum possible p-value should not be 0, but 1/N, and the max
  % should be 1
  prb_neg = prb_neg+1;
  prb_pos = prb_pos+1;
  N = Nrand+1;
end

% compute the probablities and collect the remaining details in the output structure
stat = struct();
if cfg.tail==0
  % consider both tails
  stat.prob = min(prb_neg, prb_pos)./N; % this is the probability for the most unlikely tail
  stat.stat_tfce = statobspos + statobsneg;
  stat.posdistribution = posdistribution;
  stat.negdistribution = negdistribution;
elseif cfg.tail==1
  % only consider the positive tail
  stat.prob = prb_pos./N;
  stat.stat_tfce = statobspos;
  stat.posdistribution = posdistribution;
elseif cfg.tail==-1
  % only consider the negative tail
  stat.prob = prb_neg./N;
  stat.stat_tfce = statobsneg;
  stat.negdistribution = negdistribution;
end

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

