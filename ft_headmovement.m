function [varargout] = ft_headmovement(cfg, data)

% FT_HEADMOVEMENT outputs a raw data structure, or cell-array of data structures
% reflecting the variability in the subject's head poisition relative to the 
% MEG sensors, based on continuous head position information. Current support is
% only for CTF-data. The output timeseries contain the raw HLC data, and a
% parametrization of the head movement in terms of translation and
% rotations in 3D space. The grad structure(s) have head position information
% incorporated in the coils' position/orientation and/or in the tra
% matrix, depending on the method used.
%
% Use as
%   data = ft_headmovement(cfg)
% or
%   data = ft_headmovement(cfg, data)
%
% where the configuration should contain
%   cfg.method       = string, 'updatesens' (default), 'cluster', 'avgoverrpt',
%                      'pertrial_cluster', 'pertrial' (default = 'updatesens')
% 
% and, in case no additional data structure is provided
%   cfg.dataset      = string with the filename
%
% optional arguments are
%   cfg.trl          = empty (default), or Nx3 matrix with the trial
%                      definition (see FT_DEFINETRIAL). When specified as empty,
%                      the whole recording is used. This option does not
%                      have an effect if there's an additional data input
%   cfg.numclusters  = number of segments with constant headposition in
%                      which to split the data (default = 10). This argument
%                      is only used for the methods that use clustering ('updatesens',
%                       'cluster', 'pertrial_cluster').
%
% If cfg.method = 'none', the grad in the single output structure is the grad
% as stored in the header of the recording, in dewar coordinates.
%
% If cfg.method = 'updatesens', the grad in the single output structure has
% a specification of the coils expanded as per the centroids of the position
% clusters (obtained by kmeans clustering of the HLC time series). The balancing matrix 
% is a weighted concatenation of the original tra-matrix. This method requires 
% cfg.numclusters to be specified
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
% If cfg.method = 'pertrial', the first output argument contains the time
% series of the HLC-coils, along with a grad in dewar coordinates. The
% second output argument contains a 1xNtrial cell-array of grad structures,
% expressed in head coordinates.
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
% See also FT_REGRESSCONFOUND, FT_REALTIME_HEADLOCALIZER

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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

if ft_nargin>1
  hasdata = true;
end

% set the defaults
cfg.method       = ft_getopt(cfg, 'method',      'updatesens'); % 'pertrial', 'pertrial_cluster', 'avgoverrpt', 'cluster'
cfg.numclusters  = ft_getopt(cfg, 'numclusters', 10);
cfg.feedback     = ft_getopt(cfg, 'feedback',    'yes');
cfg.computecircumcenter = ft_getopt(cfg, 'computecircumcenter', 'yes');

dokmeans = false;
if isequal(cfg.method,'updatesens') || isequal(cfg.method, 'pertrial_cluster') || isequal(cfg.method, 'cluster')
  dokmeans = true;
end

% ensure that there's data with HLC information
if ~hasdata
  % check if the input cfg is valid for this function
  cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');

  % read the header information and check whether it's a CTF dataset with HLC information
  hdr = ft_read_header(cfg.headerfile);
  assert(startsWith(ft_senstype(hdr), 'ctf'), 'currently only CTF MEG data is supported');
  assert(numel(intersect(hdr.label, {'HLC0011' 'HLC0012' 'HLC0013' 'HLC0021' 'HLC0022' 'HLC0023' 'HLC0031' 'HLC0032' 'HLC0033'}))==9, 'the data does not contain the expected head localizer channels');

  grad   = ft_datatype_sens(ctf2grad(hdr.orig, 1));
  hc     = read_ctf_hc([cfg.datafile(1:end-4),'hc']);

  % read the HLC-channels from the data file, these are expressed in m, and in dewar coordinates
  if ~isfield(cfg, 'trl') || isempty(cfg.trl)
    cfg.trl = [1 hdr.nTrials.*hdr.nSamples 0];
  end
  tmpcfg              = keepfields(cfg, {'datafile' 'trl'});
  tmpcfg.channel      = {'HLC0011' 'HLC0012' 'HLC0013' 'HLC0021' 'HLC0022' 'HLC0023' 'HLC0031' 'HLC0032' 'HLC0033'};
  tmpcfg.continuous   = 'yes';
  data                = ft_preprocessing(tmpcfg);
  data                = removefields(data, 'elec'); % this slows down a great deal
else
  % check whether the data contains HLC channel information, and whether the grad is in dewar coordinates
  sel = match_str(data.label, {'HLC0011' 'HLC0012' 'HLC0013' 'HLC0021' 'HLC0022' 'HLC0023' 'HLC0031' 'HLC0032' 'HLC0033'});
  if numel(sel)~=9
    ft_error('The required HLC channels are not present in the data');
  end
  if ~strcmp(data.grad.coordsys, 'dewar')
    ft_error('The MEG sensor array should be in dewar coordinates');
  end
  grad = data.grad;
  data = removefields(data, 'elec');
end

trial_index = cell(1,numel(data.trial));
for k = 1:numel(data.trial)
  % it sometimes happens that data are numerically 0, which causes problems downstream, replace with nans
  data.trial{k}(:,sum(data.trial{k}==0)==9) = nan;
  
  % create a bookkeeping cell-array, indexing the trial-indx
  trial_index{k} = k.*ones(1,numel(data.time{k}));
end

% average the data across time if needed by the requested method
if isequal(cfg.method, 'pertrial') || isequal(cfg.method, 'avgoverrpt') || isequal(cfg.method, 'pertrial_cluster')
  tmpcfg              = [];
  tmpcfg.avgovertime  = 'yes';
  tmpcfg.nanmean      = 'yes';
  data_timeavg        = ft_selectdata(tmpcfg, data);
  
  % concatenate across trials, each column now represents the average over time per trial
  dat  = cat(2, data_timeavg.trial{:});
  wdat = cellfun('size', data.time, 2); % weights for weighted average

  % average across trials if needed
  if isequal(cfg.method, 'avgoverrpt')
    dat = sum(dat*diag(wdat), 2)./sum(wdat);
  end
else
  % concatenate across trials, each column now represents a single sample
  dat = cat(2, data.trial{:});
end

% scale in units of the gradiometer definition, which is probably cm
dat = dat * ft_scalingfactor('m', grad.unit);

% perform the kmeans clustering if needed
if dokmeans
  if isequal(cfg.method, 'pertrial_cluster')
    trl_idx = 1:numel(data.trial);
  else
    trl_idx = cat(2, trial_index{:});

    % remove duplicates if clustering is to be performed
    [tmpdata, dum, ic] = unique(dat', 'rows');
    dat = tmpdata';

    % count how often each position occurs
    wdat = hist(ic, unique(ic));
  end

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

% find the three channels for each head localizer coil
selnas = match_str(data.label,{'HLC0011';'HLC0012';'HLC0013'});
sellpa = match_str(data.label,{'HLC0021';'HLC0022';'HLC0023'});
selrpa = match_str(data.label,{'HLC0031';'HLC0032';'HLC0033'});

if isequal(cfg.method, 'none')
  nas  = hc.dewar.nas;
  lpa  = hc.dewar.lpa;
  rpa  = hc.dewar.rpa;
else
  
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
    grad = gradnew; clear gradnew;
    grad.coordsys = 'ctf'; % hard coded
  else
    npos = size(dewar2head, 3);
    for k = 1:npos
      tmp = ft_transform_geometry(dewar2head(:,:,k), grad);
      tmp.coordsys = 'ctf';
      gradnew(k) = tmp;
    end
    grad = gradnew; clear gradnew;
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
    
  case {'avgoverrpt' 'updatesens' 'none'}
    data.grad    = grad;
    varargout{1} = data;
    
  case 'pertrial'
    varargout{1} = data;
    varargout{2} = grad;
    
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

if istrue(cfg.computecircumcenter)
  for k = 1:numel(varargout)
    for m = 1:numel(varargout{k}.trial)
      nas = varargout{k}.trial{m}(selnas,:);
      lpa = varargout{k}.trial{m}(sellpa,:);
      rpa = varargout{k}.trial{m}(selrpa,:);
      cc  = circumcenter(nas, lpa, rpa);
      varargout{k}.trial{m} = cat(1, varargout{k}.trial{m}, cc);
    end
    varargout{k}.label  = cat(1, varargout{k}.label, {'cc_xpos';'cc_ypos';'cc_zpos';'cc_xrot';'cc_yrot';'cc_zrot'});
  end
end

if istrue(cfg.feedback)
  % plot some stuff
  figure; hold on;
  title(sprintf('%s coordinates (%s)', grad.coordsys, grad.unit));
  ft_plot_axes(grad);
  ft_plot_sens(grad);
  
  fiducials = [nas;lpa;rpa];
  plot3(fiducials(:,1), fiducials(:,2), fiducials(:,3), 'b.');
  plot3(hc.dewar.nas(1), hc.dewar.nas(2), hc.dewar.nas(3), 'ro');
  plot3(hc.dewar.lpa(1), hc.dewar.lpa(2), hc.dewar.lpa(3), 'ro');
  plot3(hc.dewar.rpa(1), hc.dewar.rpa(2), hc.dewar.rpa(3), 'ro');
  axis vis3d; axis off
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble provenance
ft_postamble previous varargout
ft_postamble history varargout
ft_postamble savevar varargout

function [cc] = circumcenter(coil1, coil2, coil3)

% CIRCUMCENTER determines the position and orientation of the circumcenter
% of the three fiducial markers (MEG headposition coils).
%
% Input: X,y,z-coordinates of the 3 coils [3 X N],[3 X N],[3 X N] where N
% is timesamples/trials.
%
% Output: X,y,z-coordinates of the circumcenter [1-3 X N], and the
% orientations to the x,y,z-axes [4-6 X N].
%
% A. Stolk, 2012

% number of timesamples/trials
N = size(coil1,2);

% x-, y-, and z-coordinates of the circumcenter: use coordinates relative to point 'a' of the triangle
xba = coil2(1,:) - coil1(1,:);
yba = coil2(2,:) - coil1(2,:);
zba = coil2(3,:) - coil1(3,:);
xca = coil3(1,:) - coil1(1,:);
yca = coil3(2,:) - coil1(2,:);
zca = coil3(3,:) - coil1(3,:);

% squares of lengths of the edges incident to 'a'
balength = xba .* xba + yba .* yba + zba .* zba;
calength = xca .* xca + yca .* yca + zca .* zca;

% cross product of these edges
xcrossbc = yba .* zca - yca .* zba;
ycrossbc = zba .* xca - zca .* xba;
zcrossbc = xba .* yca - xca .* yba;

% calculate the denominator of the formulae
denominator = 0.5 ./ (xcrossbc .* xcrossbc + ycrossbc .* ycrossbc + zcrossbc .* zcrossbc);

% calculate offset (from 'a') of circumcenter
xcirca = ((balength .* yca - calength .* yba) .* zcrossbc - (balength .* zca - calength .* zba) .* ycrossbc) .* denominator;
ycirca = ((balength .* zca - calength .* zba) .* xcrossbc - (balength .* xca - calength .* xba) .* zcrossbc) .* denominator;
zcirca = ((balength .* xca - calength .* xba) .* ycrossbc - (balength .* yca - calength .* yba) .* xcrossbc) .* denominator;

% add the offset back to get the position of the origin over time.
cc(1,:) = xcirca + coil1(1,:);
cc(2,:) = ycirca + coil1(2,:);
cc(3,:) = zcirca + coil1(3,:);

% orientation of the circumcenter with respect to the x-, y-, and z-axis coordinates
v  = [cc(1,:)',    cc(2,:)',    cc(3,:)'   ];
vx = [zeros(1,N)', cc(2,:)',    cc(3,:)'   ]; % on the x-axis
vy = [cc(1,:)',    zeros(1,N)', cc(3,:)'   ]; % on the y-axis
vz = [cc(1,:)',    cc(2,:)',    zeros(1,N)']; % on the z-axis


% vectorized version of the below
normv  = sqrt(sum(v.^2,2));
thetax = acos(sum(v.*vx, 2)./(normv.*sqrt(sum(vx.^2,2))));
thetay = acos(sum(v.*vy, 2)./(normv.*sqrt(sum(vy.^2,2))));
thetaz = acos(sum(v.*vz, 2)./(normv.*sqrt(sum(vz.^2,2))));
cc(4,:) = thetax' .* (180/pi);
cc(5,:) = thetay' .* (180/pi);
cc(6,:) = thetaz' .* (180/pi);

% thetax = zeros(1,N);
% thetay = zeros(1,N);
% thetaz = zeros(1,N);
% for j = 1:N
%   % find the angles of two vectors opposing the axes
%   thetax(j) = acos(dot(v(j,:),vx(j,:))/(norm(v(j,:))*norm(vx(j,:))));
%   thetay(j) = acos(dot(v(j,:),vy(j,:))/(norm(v(j,:))*norm(vy(j,:))));
%   thetaz(j) = acos(dot(v(j,:),vz(j,:))/(norm(v(j,:))*norm(vz(j,:))));
%
%   % convert to degrees
%   cc(4,j) = (thetax(j) * (180/pi));
%   cc(5,j) = (thetay(j) * (180/pi));
%   cc(6,j) = (thetaz(j) * (180/pi));
% end

