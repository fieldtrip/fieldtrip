function [grad, grad_orig, grad_avg, trialinfo] = ft_headmovement(cfg)

% FT_HEADMOVEMENT creates a representation of the CTF gradiometer definition
% in which the headmovement in incorporated. The result can be used
% for source reconstruction.
%
% Use as
%   grad = ft_headmovement(cfg)
% where the configuration should contain
%   cfg.dataset      = string with the filename
%   cfg.trl          = Nx3 matrix with the trial definition, see FT_DEFINETRIAL
%   cfg.numclusters  = number of segments with constant headposition in which to split the data (default = 12)
%
% optional arguments are
%   cfg.updatesens   = 'yes' (default) or 'no'
%   cfg.avgovertrial = 'no' (default) or 'yes'
%   cfg.pertrial     = 'no' (default) or 'yes'
%
% If cfg.updatesens is specified, the output grad-structure has a 
% specification of the coils expanded as per the centroids of the position
% clusters. The balancing matrix is s a weighted concatenation of the
% original tra-matrix. If cfg.updatesens is false, the output grad is a
% struct-array with coil and sensor positions as per the clusters'
% centroids. If cfg.avgovertrial is specified, the clustering of head
% positions is done on the average headcoil-position per trial.
%
% Additional outputs are:
%   grad_orig = head coordinate system based grad-structure, as stored in
%                 the data directory
%   grad_avg  = head coordinate system based grad-structure, based on the 
%                 average head position as extracted from the HCL coils.
%   trialinfo = a vector ntrlx1, mapping between the grad struct-array, and
%                 the trials in cfg.trl
%
% The updatesens method and related methods are described by Stolk et al., Online and
% offline tools for head movement compensation in MEG. NeuroImage, 2012.
%
% See also FT_REGRESSCONFOUND FT_REALTIME_HEADLOCALIZER

% Copyright (C) 2008-2017, Jan-Mathijs Schoffelen, Robert Oostenveld
% Copyright (C) 2018, Jan-Mathijs Schoffelen
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

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');

% set the defaults
cfg.numclusters  = ft_getopt(cfg, 'numclusters', 12);
cfg.feedback     = ft_getopt(cfg, 'feedback',    'yes');
cfg.updatesens   = ft_getopt(cfg, 'updatesens',  'yes');
cfg.pertrial     = ft_getopt(cfg, 'pertrial',    'no');
cfg.avgovertrial = ft_getopt(cfg, 'avgovertrial', 'no');

% updatesens and pertrial are mutually exclusive
if istrue(cfg.pertrial) && istrue(cfg.updatesens)
  error('the ''pertrial'' and ''updatesens'' options are mutually exclusive');
end

% pertrial takes precendence on the clustering
if istrue(cfg.pertrial) && ~isempty(cfg.numclusters)
  cfg.numclusters = [];
end

% read the header information
hdr = ft_read_header(cfg.headerfile);
assert(numel(intersect(hdr.label, {'HLC0011' 'HLC0012' 'HLC0013' 'HLC0021' 'HLC0022' 'HLC0023' 'HLC0031' 'HLC0032' 'HLC0033'}))==9, 'the data does not contain the expected head localizer channels');

grad_head    = ctf2grad(hdr.orig, 0);
grad_head    = ft_datatype_sens(grad_head);  % ensure up-to-date sensor description (Oct 2011)
grad_dewar   = ctf2grad(hdr.orig, 1);
grad_dewar   = ft_datatype_sens(grad_dewar); % ensure up-to-date sensor description (Oct 2011)

grad = grad_dewar; % we want to work with dewar coordinates, ...
grad.chanpos = grad_head.chanpos;

% read HLC-channels
% HLC0011 HLC0012 HLC0013 x, y, z coordinates of nasion-coil in m.
% HLC0021 HLC0022 HLC0023 x, y, z coordinates of lpa-coil in m.
% HLC0031 HLC0032 HLC0033 x, y, z coordinates of rpa-coil in m.
tmpcfg              = [];
tmpcfg.dataset      = cfg.dataset;
tmpcfg.trl          = cfg.trl;
tmpcfg.channel      = {'HLC0011' 'HLC0012' 'HLC0013' 'HLC0021' 'HLC0022' 'HLC0023' 'HLC0031' 'HLC0032' 'HLC0033'};
tmpcfg.continuous   = 'yes';
data                = ft_preprocessing(tmpcfg);

if istrue(cfg.avgovertrial)
  dat  = zeros(numel(data.label), numel(data.trial));
  wdat = zeros(1, size(dat,2));
  for k = 1:numel(data.trial)
    dat(:,k)  = mean(data.trial{k},2);
    wdat(:,k) = numel(data.time{k});
  end
else
  dat = cat(2, data.trial{:});
  
  [tmpdata, ~, ic] = unique(dat', 'rows');
  dat = tmpdata';

  % count how often each position occurs
  wdat = hist(ic, unique(ic));
end
dat = dat * ft_scalingfactor('m', grad.unit); % scale in units of the gradiometer definition, which is cm

if ~isempty(cfg.numclusters)
  
  % compute the cluster means
  [bin, cluster] = kmeans(dat', cfg.numclusters);
else
  bin     = 1:size(dat,2);
  cluster = dat';
end
  
% find the three channels for each fiducial
selnas = match_str(data.label,{'HLC0011';'HLC0012';'HLC0013'});
sellpa = match_str(data.label,{'HLC0021';'HLC0022';'HLC0023'});
selrpa = match_str(data.label,{'HLC0031';'HLC0032';'HLC0033'});

ubin   = unique(bin);
nas    = zeros(numel(ubin),3);
lpa    = zeros(numel(ubin),3);
rpa    = zeros(numel(ubin),3);
numperbin = zeros(numel(ubin),1);
for k = 1:length(ubin)
  nas(k, :) = cluster(k, selnas);
  lpa(k, :) = cluster(k, sellpa);
  rpa(k, :) = cluster(k, selrpa);
  numperbin(k) = sum(wdat(bin==ubin(k)));
end
cluster_avg = sum(diag(numperbin)*cluster)./sum(numperbin);
nas_avg     = cluster_avg(:, selnas);
lpa_avg     = cluster_avg(:, sellpa);
rpa_avg     = cluster_avg(:, selrpa);


hc = read_ctf_hc([cfg.datafile(1:end-4),'hc']);
dewar2head_orig = hc.homogenous;

if istrue(cfg.feedback)
  % plot some stuff
  figure; hold on;
  title(sprintf('%s coordinates (%s)', grad_dewar.coordsys, grad_dewar.unit));
  ft_plot_axes(grad_dewar);
  ft_plot_sens(grad_dewar);
  
  fiducials = [nas;lpa;rpa];
  plot3(fiducials(:,1), fiducials(:,2), fiducials(:,3), 'b.');
  plot3(hc.dewar.nas(1), hc.dewar.nas(2), hc.dewar.nas(3), 'ro');
  plot3(hc.dewar.lpa(1), hc.dewar.lpa(2), hc.dewar.lpa(3), 'ro');
  plot3(hc.dewar.rpa(1), hc.dewar.rpa(2), hc.dewar.rpa(3), 'ro');
  axis vis3d; axis off
end

% compute transformation matrix from dewar to head coordinates
dewar2head = zeros(4, 4, size(nas,1));
for k = 1:size(dewar2head, 3)
  dewar2head(:,:,k) = ft_headcoordinates(nas(k,:), lpa(k,:), rpa(k,:), 'ctf');
end
dewar2head_avg = ft_headcoordinates(nas_avg, lpa_avg, rpa_avg, 'ctf');


if istrue(cfg.updatesens)
  npos        = size(dewar2head, 3);
  ncoils      = size(grad.coilpos,  1);
  gradnew     = grad;
  gradnew.coilpos = zeros(size(grad.coilpos,1)*npos, size(grad.coilpos,2));
  gradnew.coilori = zeros(size(grad.coilpos,1)*npos, size(grad.coilpos,2));
  gradnew.tra = repmat(grad.tra, [1 npos]);
  for m = 1:npos
    tmptransform                                  = dewar2head(:,:,m);
    gradnew.coilpos((m-1)*ncoils+1:(m*ncoils), :) = ft_warp_apply(tmptransform, grad.coilpos); % back to head coordinates
    tmptransform(1:3, 4)                          = 0; % keep rotation only
    gradnew.coilori((m-1)*ncoils+1:(m*ncoils), :) = ft_warp_apply(tmptransform, grad.coilori);
    gradnew.tra(:, (m-1)*ncoils+1:(m*ncoils))     = grad.tra.*(numperbin(m)./sum(numperbin));
  end
  
  grad = gradnew;
  trialinfo = ones(numel(data.trial),1);
else
  npos = size(dewar2head, 3);
  for k = 1:npos
    grad(k) = ft_transform_geometry(dewar2head(:,:,k), grad_dewar);
  end
  trialinfo = bin(:);
end
grad_orig = ft_transform_geometry(dewar2head_orig, grad_dewar);
grad_avg  = ft_transform_geometry(dewar2head_avg,  grad_dewar);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
