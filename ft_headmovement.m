function [varargout] = ft_headmovement(cfg)

% FT_HEADMOVEMENT creates a raw data structure, or cell-array of datastructures
% containing the HLC-coil data, which have a grad structure that has the
% head position information incorporated.
%
% Use as
%   data = ft_headmovement(cfg)
%
% where the configuration should contain
%   cfg.dataset      = string with the filename
%   cfg.method       = string, 'updatesens' (default), 'cluster', 'avgoverrpt',
%                      'pertrial_cluster', 'pertrial' (default = 'updatesens')
%
% optional arguments are
%   cfg.trl          = empty (default), or Nx3 matrix with the trial
%                      definition (see FT_DEFINETRIAL). When specified as empty,
%                      the whole recording is used.
%   cfg.numclusters  = number of segments with constant headposition in
%                      which to split the data (default = 10). This argument
%                      is only used for the clustering methods.
%
% If cfg.method = 'updatesens', the grad in the single output structure has
% a specification of the coils expanded as per the centroids of the position
% clusters. The balancing matrix is s a weighted concatenation of the
% original tra-matrix. This method requires cfg.numclusters to be specified
%
% If cfg.method = 'avgoverrpt', the grad in the single output structure has
% a specification of the coils according to the average head position
% across the specified samples.
%
% If cfg.method = 'cluster', the cell-array of output structures represent
% the epochs in which the head was considered to be positioned close to the
% corresponding kmeans-cluster's centroid. The corresponding grad-structure
% is specified according to this cluster's centroid. This method requires
% cfg.numclusters to be specified.
%
% If cfg.method = 'pertrial', the cell-array of output structures contains
% single trials, each trial with a trial-specific grad structure. Note that
% this is extremely memory inefficient with large numbers of trials, and
% probably an overkill.
%
% If cfg.method = 'pertrial_clusters', the cell-array of output structures
% contains sets of trials where the trial-specific head position was
% considered to be positioned close to the corresponding kmeans-cluster's
% centroid. The corresponding grad-structure is specified accordin to the
% cluster's centroid. This method requires cfg.numclusters to be specified.
%
% The updatesens method and related methods are described by Stolk et al., Online and
% offline tools for head movement compensation in MEG. NeuroImage, 2012.
%
% See also FT_REGRESSCONFOUND FT_REALTIME_HEADLOCALIZER

% Copyright (C) 2008-2018, Jan-Mathijs Schoffelen, Robert Oostenveld
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
% FIXME: consider allowing a data structure as input.
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');

% set the defaults
cfg.method       = ft_getopt(cfg, 'method', 'updatesens'); % 'pertrial', 'pertrial_cluster', 'avgoverrpt', 'cluster'
cfg.numclusters  = ft_getopt(cfg, 'numclusters', 10);
cfg.feedback     = ft_getopt(cfg, 'feedback',    'yes');

if isequal(cfg.method,'updatesens') || isequal(cfg.method, 'pertrial_cluster') || isequal(cfg.method, 'cluster')
  dokmeans = true;
else
  dokmeans = false;
end

% read the header information and check whether it's a CTF dataset with HLC information
hdr = ft_read_header(cfg.headerfile);
assert(numel(intersect(hdr.label, {'HLC0011' 'HLC0012' 'HLC0013' 'HLC0021' 'HLC0022' 'HLC0023' 'HLC0031' 'HLC0032' 'HLC0033'}))==9, 'the data does not contain the expected head localizer channels');

grad_head    = ctf2grad(hdr.orig, 0);
grad_head    = ft_datatype_sens(grad_head);  % ensure up-to-date sensor description (Oct 2011)
grad_dewar   = ctf2grad(hdr.orig, 1);
grad_dewar   = ft_datatype_sens(grad_dewar); % ensure up-to-date sensor description (Oct 2011)

grad = grad_dewar; % we want to work with dewar coordinates, ...
grad.chanpos = grad_head.chanpos;

% read the HLC-channels
% HLC0011 HLC0012 HLC0013 x, y, z coordinates of nasion-coil in m.
% HLC0021 HLC0022 HLC0023 x, y, z coordinates of lpa-coil in m.
% HLC0031 HLC0032 HLC0033 x, y, z coordinates of rpa-coil in m.
if ~isfield(cfg, 'trl') || isempty(cfg.trl)
  cfg.trl = [1 hdr.nTrials.*hdr.nSamples 0];
end
tmpcfg              = [];
tmpcfg.dataset      = cfg.dataset;
tmpcfg.trl          = cfg.trl;
tmpcfg.channel      = {'HLC0011' 'HLC0012' 'HLC0013' 'HLC0021' 'HLC0022' 'HLC0023' 'HLC0031' 'HLC0032' 'HLC0033'};
tmpcfg.continuous   = 'yes';
data                = ft_preprocessing(tmpcfg);
data                = removefields(data, 'elec'); % this slows down a great
% rendering the persistent variable trick useless.
% we don't need the elec anyway

wdat                = cellfun('size', data.time, 2); % weights for weighted average

trial_index = cell(1,numel(data.trial));
for k = 1:numel(data.trial)
  % it sometimes happens that data are numerically 0, which causes problems downstream, replace with nans
  data.trial{k}(:,sum(data.trial{k}==0)==9) = nan;
  
  % create a bookkeeping cell-array, indexing the trial-indx
  trial_index{k} = k.*ones(1,numel(data.time{k}));
end

% average across time if needed
if isequal(cfg.method, 'pertrial') || isequal(cfg.method, 'avgoverrpt') || isequal(cfg.method, 'pertrial_cluster')
  tmpcfg              = [];
  tmpcfg.avgovertime  = 'yes';
  tmpcfg.nanmean      = 'yes';
  data_timeavg        = ft_selectdata(tmpcfg, data);
  
  % concatenate across trials and scale the units
  dat = cat(2, data_timeavg.trial{:});
else
  % concatenate across trials and scale the units
  dat = cat(2, data.trial{:});
end

% scale in units of the gradiometer definition, which is probably cm
dat = dat * ft_scalingfactor('m', grad.unit);

if isequal(cfg.method, 'pertrial_cluster')
  trl_idx = 1:numel(data.trial);
else
  trl_idx = cat(2, trial_index{:});
end

% average across trials if needed
if isequal(cfg.method, 'avgoverrpt')
  dat = sum(dat*diag(wdat), 2)./sum(wdat);
end

% remove duplicates if clustering is to be performed
if dokmeans && ~isequal(cfg.method, 'pertrial_cluster')
  [tmpdata, dum, ic] = unique(dat', 'rows');
  dat = tmpdata';
  
  % count how often each position occurs
  wdat = hist(ic, unique(ic));
end

% perform the clustering if needed
if dokmeans
  % compute the cluster means
  [bin, dat] = kmeans(dat', cfg.numclusters, 'EmptyAction', 'drop');
  
  % create a cell-array 1xnrpt with time specific indices of cluster id
  cluster_id = cell(1,numel(data.trial));
  for k = 1:numel(data.trial)
    cluster_id{k} = nan+zeros(1,numel(data.time{k}));
    if ~isequal(cfg.method, 'pertrial_cluster')
      for m = 1:size(dat,1)
        tmpdat = ic(trl_idx==k);
        cluster_id{k}(ismember(tmpdat, find(bin==m))) = m;
      end
    else
      cluster_id{k}(:) = bin(k);
    end
  end
  
else
  bin = 1:size(dat,2);
  dat = dat';
end

% find the three channels for each fiducial
selnas = match_str(data.label,{'HLC0011';'HLC0012';'HLC0013'});
sellpa = match_str(data.label,{'HLC0021';'HLC0022';'HLC0023'});
selrpa = match_str(data.label,{'HLC0031';'HLC0032';'HLC0033'});

ubin   = unique(bin(isfinite(bin)));
nas    = zeros(numel(ubin),3);
lpa    = zeros(numel(ubin),3);
rpa    = zeros(numel(ubin),3);
numperbin = zeros(numel(ubin),1);
for k = 1:length(ubin)
  nas(k, :) = dat(k, selnas);
  lpa(k, :) = dat(k, sellpa);
  rpa(k, :) = dat(k, selrpa);
  numperbin(k) = sum(wdat(bin==ubin(k)));
end

hc = read_ctf_hc([cfg.datafile(1:end-4),'hc']);

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

if isequal(cfg.method, 'updatesens')
  npos        = size(dewar2head, 3);
  ncoils      = size(grad.coilpos,  1);
  gradnew     = grad;
  gradnew.coilpos = zeros(size(grad.coilpos,1)*npos, size(grad.coilpos,2));
  gradnew.coilori = zeros(size(grad.coilpos,1)*npos, size(grad.coilpos,2));
  gradnew.tra = repmat(grad.tra, [1 npos]);
  for m = 1:npos
    tmptransform                                  = dewar2head(:,:,m);
    gradnew.coilpos((m-1)*ncoils+1:(m*ncoils), :) = ft_warp_apply(tmptransform, grad.coilpos); % back to head coordinates
    tmptransform(1:3, 4)                          = 0; % keep only the rotation
    gradnew.coilori((m-1)*ncoils+1:(m*ncoils), :) = ft_warp_apply(tmptransform, grad.coilori);
    gradnew.tra(:, (m-1)*ncoils+1:(m*ncoils))     = grad.tra.*(numperbin(m)./sum(numperbin));
  end
  grad = gradnew;
else
  npos = size(dewar2head, 3);
  for k = 1:npos
    grad(k) = ft_transform_geometry(dewar2head(:,:,k), grad_dewar);
  end
end

% prepare the output data
switch cfg.method
  case 'cluster'
    varargout     = cell(1,numel(grad));
    tmpdata       = data;
    tmpdata.trial = cluster_id;
    tmpdata.label = {'cluster_id'};
    data          = ft_appenddata([],data,tmpdata);
    
    for k = 1:numel(grad)
      tmpcfg = [];
      tmpcfg.artfctdef.bpfilter = 'no';
      tmpcfg.artfctdef.threshold.channel   = {'cluster_id'};
      tmpcfg.artfctdef.threshold.min       = 0.9+k-1;
      tmpcfg.artfctdef.threshold.max       = 1.1+k-1;
      tmpcfg.artfctdef.threshold.bpfilter  = 'no';
      tmpcfg = ft_artifact_threshold(tmpcfg, tmpdata);
      artifacts = tmpcfg.artfctdef.threshold.artifact;
      
      tmpcfg = [];
      tmpcfg.artfctdef.reject = 'partial';
      tmpcfg.artfctdef.threshold.artifact = artifacts;
      tmpdata_clus = ft_rejectartifact(tmpcfg, data);
      tmpdata_clus.grad = grad(k);
      
      varargout{k} = tmpdata_clus;
    end
    
  case {'avgoverrpt' 'updatesens'}
    data.grad    = grad;
    varargout{1} = data;
    
  case 'pertrial'
    ft_error('not yet implemented');
    
  case 'pertrial_cluster'
    varargout     = cell(1,numel(grad));
    tmpdata       = data;
    tmpdata.trial = cluster_id;
    tmpdata.label = {'cluster_id'};
    data          = ft_appenddata([],data,tmpdata);
    for k = 1:numel(grad)
      tmpcfg        = [];
      tmpcfg.trials = find(bin==k);
      tmpdata_clus  = ft_selectdata(tmpcfg, data);
      %tmpcfg.previous = tmpdata_clus.cfg;
      
      %tmpdata_clus.cfg  = tmpcfg;
      tmpdata_clus.grad = grad(k);
      
      varargout{k} = tmpdata_clus;
    end
end % switch method

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous varargout
ft_postamble history varargout
ft_postamble savevar varargout

