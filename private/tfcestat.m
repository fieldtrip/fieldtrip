function [stattfce, cfg] = tfcestat(cfg, statobs)

% TFCESTAT computes threshold-free cluster statistic for multidimensional
% channel-freq-time or volumetric source data.
%
% Two implementations are available, selected with cfg.tfce_method:
%
%   'exact'    (default) evaluates the TFCE integral EXACTLY (eTFCE). Because
%              the cluster extent e(h) is piecewise constant in the height h
%              (it changes only when h crosses a data value), the integral
%
%                  TFCE(v) = integral_{h0}^{h_v} e(h)^E * h^H dh
%
%              has a closed form on each interval and is obtained in a single
%              pass with a disjoint-set (union-find) cluster-retrieval forest,
%              using the same spatial connectivity as FINDCLUSTER. This removes
%              the cfg.tfce_nsteps discretisation (and its bias), is typically
%              much faster than the threshold loop, and needs no SPM toolbox.
%
%   'discrete' the original implementation, which approximates the integral by
%              summing over cfg.tfce_nsteps discrete height thresholds and
%              calls FINDCLUSTER (requires SPM) at each threshold. Kept for
%              backward compatibility and for validating the exact method.
%
% Relevant options:
%   cfg.tfce_method = 'exact' (default) or 'discrete'
%   cfg.tfce_H      = height exponent H (default 2)
%   cfg.tfce_E      = extent exponent E (default 0.5)
%   cfg.tfce_h0     = baseline/offset that is subtracted first (default 0)
%   cfg.tfce_nsteps = number of height steps for the 'discrete' method (default 100)
%   cfg.tail        = -1, 0 or 1 (default 0)
%
% The exact-TFCE (eTFCE) algorithm was developed by Xu Chen, Wouter D. Weeda,
% Thomas E. Nichols and Jelle J. Goeman (2026), "eTFCE: Exact Threshold-Free
% Cluster Enhancement via Fast Cluster Retrieval", arXiv:2603.03004. The 'exact'
% method below is a MATLAB implementation of that published algorithm for
% FieldTrip; it was NOT created or developed by the implementer, who only
% translated the method into code.
%
% See also CLUSTERSTAT, FINDCLUSTER

% Copyright (C) 2021, Jan-Mathijs Schoffelen (discretised implementation)
%
% The exact (eTFCE) algorithm was created by:
% Copyright (C) 2026, Xu Chen, Wouter D. Weeda, Thomas E. Nichols and Jelle J. Goeman (eTFCE algorithm, arXiv:2603.03004)
% and implemented for FieldTrip (the implementer did not develop the method) by:
% Copyright (C) 2026, Devon Yanitski (FieldTrip implementation of the eTFCE algorithm)
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
% honour the documented cfg.tfce_nsteps (set by the caller), falling back to the
% legacy 'nsteps' field and finally to 100; only used by the 'discrete' method
cfg.tfce_nsteps  = ft_getopt(cfg, 'tfce_nsteps',  ft_getopt(cfg, 'nsteps', 100));
cfg.tfce_method  = ft_getopt(cfg, 'tfce_method',  'exact');   % 'exact' (eTFCE) or 'discrete'
cfg.height       = ft_getopt(cfg, 'height',       []);
% these defaults are already set in the caller function,
% but may be necessary if a user calls this function directly
cfg.connectivity = ft_getopt(cfg, 'connectivity', false);

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

% compute the height (max statistic) once and return it in cfg, so that in a
% subsequent call (for the randomizations) the same value is reused. The
% 'discrete' method needs it for the stepsize; the 'exact' method does not need
% it for the integral, but it keeps the cfg round-trip identical.
if isempty(cfg.height)
  if needneg && needpos
    cfg.height = max(abs(statobs(:)));
  elseif needneg
    cfg.height = max(-statobs(:));
  elseif needpos
    cfg.height = max(statobs(:));
  end
end

switch lower(cfg.tfce_method)
  case 'exact'
    stattfce = tfce_exact(cfg, statobs, connmat, needpos, needneg);
  case 'discrete'
    % ensure that the preferred SPM version is on the path (findcluster needs it)
    ft_hastoolbox(cfg.spmversion, 1);
    stattfce = tfce_discrete(cfg, statobs, connmat, needpos, needneg);
  otherwise
    ft_error('unsupported cfg.tfce_method ''%s'' (use ''exact'' or ''discrete'')', cfg.tfce_method);
end


%==========================================================================
% DISCRETE method: original implementation, summing over cfg.tfce_nsteps
% height thresholds and calling findcluster at each one.
%==========================================================================
function stattfce = tfce_discrete(cfg, statobs, connmat, needpos, needneg)

stepsize = cfg.height./cfg.tfce_nsteps;

% first do the clustering on the observed data
spacereshapeable = (isscalar(connmat) && ~isfinite(connmat));

statobspos = 0;
statobsneg = 0;

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


%==========================================================================
% EXACT method (eTFCE): evaluate the TFCE integral in a single pass with a
% union-find cluster-retrieval forest, using the same connectivity as the
% discrete method but without thresholding at discrete heights.
%==========================================================================
function stattfce = tfce_exact(cfg, statobs, connmat, needpos, needneg)

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
