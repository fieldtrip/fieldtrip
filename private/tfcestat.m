function [stattfce, cfg] = tfcestat(cfg, statobs)

% TFCESTAT computes threshold-free cluster statistic for multidimensional
% channel-freq-time or volumetric source data.
%
% This implementation evaluates the TFCE integral EXACTLY (eTFCE) rather than
% approximating it with a fixed number of discrete height thresholds. Because
% the cluster extent e(h) is piecewise constant in the height h (it changes
% only when h crosses a data value), the integral
%
%     TFCE(v) = integral_{h0}^{h_v} e(h)^E * h^H dh
%
% has a closed form on each interval and is obtained in a single pass with a
% disjoint-set (union-find) cluster-retrieval forest, using the same spatial
% connectivity as FINDCLUSTER. This removes the cfg.tfce_nsteps discretisation
% (and its bias), is typically much faster than the threshold loop, and needs
% no SPM/Image-Processing clustering routine.
%
% The exact-TFCE (eTFCE) algorithm is the method of Chen, Weeda, Nichols &
% Goeman (2026), "eTFCE: Exact Threshold-Free Cluster Enhancement via Fast
% Cluster Retrieval", arXiv:2603.03004; this is a MATLAB implementation of it
% for FieldTrip and was not developed by the FieldTrip team.
%
% See also CLUSTERSTAT, FINDCLUSTER

% Copyright (C) 2021, Jan-Mathijs Schoffelen (original discretised version)
% Copyright (C) 2026, exact (eTFCE) reimplementation
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
cfg.tfce_nsteps  = ft_getopt(cfg, 'nsteps',       100);       % unused by the exact method; kept for cfg compatibility
cfg.height       = ft_getopt(cfg, 'height',       []);
% these defaults are already set in the caller function,
% but may be necessary if a user calls this function directly
cfg.connectivity = ft_getopt(cfg, 'connectivity', false);

% NB the exact method does not threshold at a discrete set of heights, so it
% does not call findcluster/spm_bwlabel and therefore needs no SPM toolbox.

if isempty(cfg.dim)
  ft_error('cfg.dim should be defined and not empty');
end

if isempty(cfg.inside)
  cfg.inside = true(cfg.dim);
end % cfg.inside is set in ft_sourcestatistics, but is also needed for timelock and freq

if isfield(cfg, 'origdim')
  cfg.dim = cfg.origdim;
end % this snippet is to support correct clustering of N-dimensional data

% get connectivity matrix for the spatially neighbouring elements
connmat = full(ft_getopt(cfg, 'connectivity', false));

needpos = cfg.tail==0 || cfg.tail== 1;
needneg = cfg.tail==0 || cfg.tail==-1;

% remove the offset, which by default is 0
statobs = statobs - cfg.tfce_h0;

% explicitly compute the height (max statistic), to be returned in cfg so that
% subsequent (randomization) calls can reuse the same value; the exact method
% does not need it for the integral but it keeps the cfg round-trip unchanged.
if isempty(cfg.height)
  if needneg && needpos
    cfg.height = max(abs(statobs(:)));
  elseif needneg
    cfg.height = max(-statobs(:));
  elseif needpos
    cfg.height = max(statobs(:));
  end
end

% layout: replicate how the discretised code arranges the data for findcluster
spacereshapeable = (isscalar(connmat) && ~isfinite(connmat));
if spacereshapeable
  % source data on a 3D grid: a fake leading singleton spatial dimension
  arrsz = [1 cfg.dim];
  fullvals = zeros(prod(arrsz),1);
  fullvals(cfg.inside(:)) = statobs(:);
else
  % channel x (freq) x time data: spatial (channel) dimension is first
  arrsz = [cfg.dim 1];
  fullvals = statobs(:);
end
Nlin = prod(arrsz);

% build the voxel adjacency once (geometry only), matching findcluster:
% face connectivity over the non-spatial dimensions + connmat across the
% spatial (first) dimension at matching non-spatial positions
edges = local_build_edges(arrsz, connmat);

% exact positive and negative TFCE
statobspos = zeros(numel(statobs),1);
statobsneg = zeros(numel(statobs),1);

if needpos
  sc = local_etfce(max(fullvals,0), edges, Nlin, cfg.tfce_E, cfg.tfce_H);
  if spacereshapeable, statobspos = sc(cfg.inside(:)); else, statobspos = sc(:); end
end

if needneg
  sc = local_etfce(max(-fullvals,0), edges, Nlin, cfg.tfce_E, cfg.tfce_H);
  if spacereshapeable, statobsneg = -sc(cfg.inside(:)); else, statobsneg = -sc(:); end
end

stattfce = statobspos + statobsneg;


%==========================================================================
% exact TFCE for the (non-negative) input vals, returned per voxel
%==========================================================================
function score = local_etfce(vals, edges, Nlin, E, H)

score  = zeros(Nlin,1);
active = vals > 0;
N      = nnz(active);
if N==0, return; end

activeIdx         = find(active);
nodeid            = zeros(Nlin,1);
nodeid(activeIdx) = 1:N;
nodeval           = vals(activeIdx);

ea = nodeid(edges(:,1));
eb = nodeid(edges(:,2));
keep = ea>0 & eb>0;
ea = ea(keep); eb = eb(keep);

src = [ea; eb]; dst = [eb; ea];
[src, o] = sort(src); dst = dst(o);
cnt    = accumarray(src, 1, [N 1]);
startp = cumsum([1; cnt]);

[~, order] = sort(nodeval, 'descend');
rnk        = zeros(N,1); rnk(order) = 1:N;

ufp = (1:N)'; csize = ones(N,1); croot = (1:N)'; hpar = (1:N)'; ext = ones(N,1);
for p = 1:N
  i = order(p); s = startp(i); e = startp(i+1)-1;
  ri = i; while ufp(ri)~=ri, ri = ufp(ri); end
  a = i; while ufp(a)~=ri, t = ufp(a); ufp(a) = ri; a = t; end
  for k = s:e
    j = dst(k);
    if rnk(j) >= p, continue; end
    rj = j; while ufp(rj)~=rj, rj = ufp(rj); end
    a = j; while ufp(a)~=rj, t = ufp(a); ufp(a) = rj; a = t; end
    if rj==ri, continue; end
    cr = croot(rj); hpar(cr) = i;
    if csize(ri) >= csize(rj)
      ufp(rj) = ri; csize(ri) = csize(ri)+csize(rj); croot(ri) = i;
    else
      ufp(ri) = rj; csize(rj) = csize(ri)+csize(rj); croot(rj) = i; ri = rj;
    end
  end
  ext(i) = csize(ri);
end

T = zeros(N,1); c = 1/(H+1);
for p = N:-1:1
  i = order(p); u = hpar(i);
  if u==i
    T(i) = c*(ext(i)^E)*(nodeval(i)^(H+1));
  else
    T(i) = T(u) + c*(ext(i)^E)*(nodeval(i)^(H+1) - nodeval(u)^(H+1));
  end
end
score(activeIdx) = T;


%==========================================================================
% voxel adjacency matching FINDCLUSTER, as an [M x 2] edge list of linear
% indices into an array of size arrsz (spatial/channel dimension first)
%==========================================================================
function edges = local_build_edges(arrsz, connmat)

D1    = arrsz(1);
nd    = numel(arrsz);
ngrid = prod(arrsz(2:end));
edges = zeros(0,2);

% face-connected grid edges over the non-spatial dimensions
A = reshape(1:prod(arrsz), arrsz);
for d = 2:nd
  if arrsz(d) >= 2
    S1 = repmat({':'},1,nd); S1{d} = 1:arrsz(d)-1;
    S2 = repmat({':'},1,nd); S2{d} = 2:arrsz(d);
    a = A(S1{:}); b = A(S2{:});
    edges = [edges; a(:), b(:)]; %#ok<AGROW>
  end
end

% spatial (channel) edges from connmat, at matching non-spatial positions
if ~isscalar(connmat) && all(size(connmat)==D1)
  M = (connmat~=0); M = triu(M | M', 1);
  [c1,c2] = find(M);
  if ~isempty(c1)
    poff = (0:ngrid-1)*D1;            % 1 x ngrid offsets
    a = c1(:) + poff;                 % P x ngrid
    b = c2(:) + poff;
    edges = [edges; a(:), b(:)];
  end
end
