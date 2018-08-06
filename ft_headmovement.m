function [grad] = ft_headmovement(cfg)

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
% This method and related methods are described by Stolk et al., Online and
% offline tools for head movement compensation in MEG. NeuroImage, 2012.
%
% See also FT_REGRESSCONFOUND FT_REALTIME_HEADLOCALIZER

% Copyright (C) 2008-2017, Jan-Mathijs Schoffelen, Robert Oostenveld
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
cfg.numclusters = ft_getopt(cfg, 'numclusters', 12);
cfg.feedback    = ft_getopt(cfg, 'feedback', 'yes');

% read the header information
hdr = ft_read_header(cfg.headerfile);
assert(numel(intersect(hdr.label, {'HLC0011' 'HLC0012' 'HLC0013' 'HLC0021' 'HLC0022' 'HLC0023' 'HLC0031' 'HLC0032' 'HLC0033'}))==9, 'the data does not contain the expected head localizer channels');

% work with gradiometers in dewar coordinates, since HLCs are also
% in dewar coords. at present I did not find the nas, lpa, rpa channels,
% which according to ctf's documentation should contain the positions
% of these channels directly (HDAC channels). FIXME
grad_head    = ctf2grad(hdr.orig, 0);
grad_dewar   = ctf2grad(hdr.orig, 1);
grad_head    = ft_datatype_sens(grad_head);  % ensure up-to-date sensor description (Oct 2011)
grad_dewar   = ft_datatype_sens(grad_dewar); % ensure up-to-date sensor description (Oct 2011)

grad         = grad_dewar;        % we want to work with dewar coordinates, ...
grad.chanpos = grad_head.chanpos; % except the chanpos, which should remain in head coordinates

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

dat = cat(2, data.trial{:});
dat = dat * ft_scalingfactor('m', grad.unit); % scale in units of the gradiometer definition

[tmpdata, dum, ic] = unique(dat', 'rows');
dat = tmpdata';

% counthow often each position occurs
[wdat, dum] = hist(ic, unique(ic));

% compute the cluster means
[bin, cluster] = kmeans(dat', cfg.numclusters);

% find the three channels for each fiducial
selnas = strmatch('HLC001', data.label);
sellpa = strmatch('HLC002', data.label);
selrpa = strmatch('HLC003', data.label);

ubin   = unique(bin);
for k = 1:length(ubin)
  nas(k, :) = cluster(k, selnas);
  lpa(k, :) = cluster(k, sellpa);
  rpa(k, :) = cluster(k, selrpa);
  numperbin(k) = sum(wdat(bin==ubin(k)));
end

if istrue(cfg.feedback)
  hc    = read_ctf_hc([cfg.datafile(1:end-4),'hc']);
  
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
transform = zeros(4, 4, size(nas,1));
for k = 1:size(transform, 3)
  transform(:,:,k) = ft_headcoordinates(nas(k,:), lpa(k,:), rpa(k,:), 'ctf');
end

npos        = size(transform, 3);
ncoils      = size(grad.coilpos,  1);
gradnew     = grad;
gradnew.coilpos = zeros(size(grad.coilpos,1)*npos, size(grad.coilpos,2));
gradnew.coilori = zeros(size(grad.coilpos,1)*npos, size(grad.coilpos,2));
gradnew.tra = repmat(grad.tra, [1 npos]);
for m = 1:npos
  tmptransform                                  = transform(:,:,m);
  gradnew.coilpos((m-1)*ncoils+1:(m*ncoils), :) = ft_warp_apply(tmptransform, grad.coilpos); % back to head coordinates
  tmptransform(1:3, 4)                          = 0; % keep rotation only
  gradnew.coilori((m-1)*ncoils+1:(m*ncoils), :) = ft_warp_apply(tmptransform, grad.coilori);
  gradnew.tra(:, (m-1)*ncoils+1:(m*ncoils))     = grad.tra.*(numperbin(m)./sum(numperbin));
end

grad = gradnew;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
