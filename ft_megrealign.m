function [data] = ft_megrealign(cfg, data)

% FT_MEGREALIGN interpolates MEG data towards standard gradiometer locations by
% projecting the individual timelocked data towards a coarse source reconstructed
% representation and computing the magnetic field on the standard gradiometer
% locations.
%
% Use as
%   [interp] = ft_megrealign(cfg, data)
%
% Required configuration options:
%   cfg.template
%   cfg.inwardshift
%
% The new gradiometer definition is obtained from a template dataset,
% or can be constructed by averaging the gradiometer positions over
% multiple datasets.
%   cfg.template       = single dataset that serves as template
%   cfg.template(1..N) = datasets that are averaged into the standard
%
% The realignment is done by computing a minumum norm estimate using a
% large number of dipoles that are placed in the upper layer of the brain
% surface, followed by a forward computation towards the template
% gradiometer array. This requires the specification of a volume conduction
% model of the head and of a source model.
%
% A volume conduction model of the head should be specified with
%   cfg.headmodel   = structure, see FT_PREPARE_HEADMODEL
%
% A source model (i.e. a superficial layer with distributed sources) can be
% constructed from a headshape file, or from inner surface of the volume conduction
% model using FT_PREPARE_SOIURCEMODEL using the following options
%   cfg.spheremesh  = number of dipoles in the source layer (default = 642)
%   cfg.inwardshift = depth of the source layer relative to the headshape
%                     surface or volume conduction model (no default
%                     supplied, see below)
%   cfg.headshape   = a filename containing headshape, a structure containing a
%                     single triangulated boundary, or a Nx3 matrix with surface
%                     points
%
% If you specify a headshape and it describes the skin surface, you should specify an
% inward shift of 2.5 cm.
%
% For a single-sphere or a local-spheres volume conduction model based on the skin
% surface, an inward shift of 2.5 cm is reasonable.
%
% For a single-sphere or a local-spheres volume conduction model based on the brain
% surface, you should probably use an inward shift of about 1 cm.
%
% For a realistic single-shell volume conduction model based on the brain surface, you
% should probably use an inward shift of about 1 cm.
%
% Other options are
% cfg.pruneratio  = for singular values, default is 1e-3
% cfg.verify      = 'yes' or 'no', show the percentage difference (default = 'yes')
% cfg.feedback    = 'yes' or 'no' (default = 'no')
% cfg.channel     =  Nx1 cell-array with selection of channels (default = 'MEG'),
%                      see FT_CHANNELSELECTION for details
% cfg.trials      = 'all' or a selection given as a 1xN vector (default = 'all')
%
% This implements the method described by T.R. Knosche, Transformation
% of whole-head MEG recordings between different sensor positions.
% Biomed Tech (Berl). 2002 Mar;47(3):59-62. For more information and
% related methods, see Stolk et al., Online and offline tools for head
% movement compensation in MEG. NeuroImage, 2012.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_PREPARE_LOCALSPHERES, FT_PREPARE_SINGLESHELL

% Copyright (C) 2004-2014, Robert Oostenveld
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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',    {'plot3d',      'feedback'});
cfg = ft_checkconfig(cfg, 'renamedval', {'headshape',   'headmodel', []});
cfg = ft_checkconfig(cfg, 'required',   {'inwardshift', 'template'});
cfg = ft_checkconfig(cfg, 'renamed',    {'hdmfile',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'vol',         'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'grid',        'sourcemodel'});

% set the default configuration
cfg.headshape  = ft_getopt(cfg, 'headshape', []);
cfg.pruneratio = ft_getopt(cfg, 'pruneratio', 1e-3);
cfg.spheremesh = ft_getopt(cfg, 'spheremesh', 642);
cfg.verify     = ft_getopt(cfg, 'verify',     'yes');
cfg.feedback   = ft_getopt(cfg, 'feedback',   'yes');
cfg.trials     = ft_getopt(cfg, 'trials',     'all', 1);
cfg.channel    = ft_getopt(cfg, 'channel',    'MEG');
cfg.topoparam  = ft_getopt(cfg, 'topoparam',  'rms');

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes', 'ismeg', 'yes');

% do realignment per trial
pertrial = all(ismember({'nasX';'nasY';'nasZ';'lpaX';'lpaY';'lpaZ';'rpaX';'rpaY';'rpaZ'}, data.label));

% put the low-level options pertaining to the dipole grid in their own field
cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'}); % this is moved to cfg.sourcemodel.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.sourcemodel.unit by the subsequent createsubcfg

% put the low-level options pertaining to the sourcemodel in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg', {'sourcemodel'});
% move some fields from cfg.sourcemodel back to the top-level configuration
cfg = ft_checkconfig(cfg, 'createtopcfg', {'sourcemodel'});

if isstruct(cfg.template)
  % this should be a cell-array
  cfg.template = {cfg.template};
end

% retain only the MEG channels in the data and temporarily store
% the rest, these will be added back to the transformed data later.

% select trials and channels of interest
tmpcfg = [];
tmpcfg.trials  = cfg.trials;
tmpcfg.channel = setdiff(data.label, ft_channelselection(cfg.channel, data.label));
rest = ft_selectdata(tmpcfg, data);

tmpcfg.channel = ft_channelselection(cfg.channel, data.label);
data = ft_selectdata(tmpcfg, data);

% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

Ntrials = length(data.trial);

% cfg.channel = ft_channelselection(cfg.channel, data.label);
% dataindx = match_str(data.label, cfg.channel);
% restindx = setdiff(1:length(data.label),dataindx);
% if ~isempty(restindx)
%   fprintf('removing %d non-MEG channels from the data\n', length(restindx));
%   rest.label = data.label(restindx);    % first remember the rest
%   data.label = data.label(dataindx);    % then reduce the data
%   for i=1:Ntrials
%     rest.trial{i} = data.trial{i}(restindx,:);  % first remember the rest
%     data.trial{i} = data.trial{i}(dataindx,:);  % then reduce the data
%   end
% else
%   rest.label = {};
%   rest.trial = {};
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the average template gradiometer array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
template = struct([]); % initialize as empty structure
for i=1:length(cfg.template)
  if ischar(cfg.template{i})
    fprintf('reading template sensor position from %s\n', cfg.template{i});
    tmp = ft_read_sens(cfg.template{i}, 'senstype', 'meg');
  elseif isstruct(cfg.template{i}) && isfield(cfg.template{i}, 'coilpos') && isfield(cfg.template{i}, 'coilori') && isfield(cfg.template{i}, 'tra')
    tmp = cfg.template{i};
  elseif isstruct(cfg.template{i}) && isfield(cfg.template{i}, 'pnt') && isfield(cfg.template{i}, 'ori') && isfield(cfg.template{i}, 'tra')
    % it seems to be a pre-2011v1 type gradiometer structure, update it
    tmp = ft_datatype_sens(cfg.template{i});
  else
    ft_error('unrecognized template input');
  end
  % prevent "Subscripted assignment between dissimilar structures" error
  template = appendstruct(template, tmp); clear tmp
end

grad = ft_average_sens(template);

% construct the final template gradiometer definition
template = [];
template.grad = grad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FT_PREPARE_VOL_SENS will match the data labels, the gradiometer labels and the
% volume model labels (in case of a localspheres model) and result in a gradiometer
% definition that only contains the gradiometers that are present in the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volcfg = [];
volcfg.headmodel = cfg.headmodel;
volcfg.grad      = data.grad;
volcfg.channel   = data.label; % this might be a subset of the MEG channels
gradorig         = data.grad;  % this is needed later on for plotting. As of
% yet the next step is not entirely correct, because it does not keep track
% of the balancing of the gradiometer array. FIXME this may require some
% thought because the leadfields are computed with low level functions and
% do not easily accommodate for matching the correct channels with each
% other (in order to compute the projection matrix).
[volold, data.grad] = prepare_headmodel(volcfg);

% note that it is neccessary to keep the two volume conduction models
% separate, since the single-shell Nolte model contains gradiometer specific
% precomputed parameters. Note that this is not guaranteed to result in a
% good projection for local sphere models.
volcfg.grad    = template.grad;
volcfg.channel = 'MEG'; % include all MEG channels
[volnew, template.grad] = prepare_headmodel(volcfg);

if strcmp(ft_senstype(data.grad), ft_senstype(template.grad))
  [id, it] = match_str(data.grad.label, template.grad.label);
  fprintf('mean distance towards template gradiometers is %.2f %s\n', mean(sum((data.grad.chanpos(id,:)-template.grad.chanpos(it,:)).^2, 2).^0.5), template.grad.unit);
else
  % the projection is from one MEG system to another MEG system, which makes a comparison of the data difficult
  cfg.feedback = 'no';
  cfg.verify = 'no';
end

% copy all options that are potentially used in ft_prepare_sourcemodel
tmpcfg           = keepfields(cfg, {'sourcemodel', 'mri', 'headshape', 'symmetry', 'smooth', 'threshold', 'spheremesh', 'inwardshift', 'xgrid' 'ygrid', 'zgrid', 'resolution', 'tight', 'warpmni', 'template', 'showcallinfo'});
tmpcfg.headmodel = volold;
tmpcfg.grad      = data.grad;
% create the source positions on which the data will be projected
sourcemodel = ft_prepare_sourcemodel(tmpcfg);
pos = sourcemodel.pos;

% sometimes some of the dipole positions are nan, due to problems with the headsurface triangulation
% remove them to prevent problems with the forward computation
sel = find(any(isnan(pos(:,1)),2));
pos(sel,:) = [];

% compute the forward model for the new gradiometer positions
fprintf('computing forward model for %d dipoles\n', size(pos,1));
lfnew = ft_compute_leadfield(pos, template.grad, volnew);
if ~pertrial
  %this needs to be done only once
  lfold = ft_compute_leadfield(pos, data.grad, volold);
  [realign, noalign, bkalign] = computeprojection(lfold, lfnew, cfg.pruneratio, cfg.verify);
else
  %the forward model and realignment matrices have to be computed for each trial
  %this also goes for the singleshell volume conductor model
  %x = which('rigidbodyJM'); %this function is needed
  %if isempty(x)
  %  ft_error('you are trying out experimental code for which you need some extra functionality which is currently not in the release version of FieldTrip. if you are interested in trying it out, contact Jan-Mathijs');
  %end
end

% interpolate the data towards the template gradiometers
for i=1:Ntrials
  fprintf('realigning trial %d\n', i);
  if pertrial
    %warp the gradiometer array according to the motiontracking data
    sel   = match_str(rest.label, {'nasX';'nasY';'nasZ';'lpaX';'lpaY';'lpaZ';'rpaX';'rpaY';'rpaZ'});
    hmdat = rest.trial{i}(sel,:);
    if ~all(hmdat==repmat(hmdat(:,1),[1 size(hmdat,2)]))
      ft_error('only one position per trial is at present allowed');
    else
      %M    = rigidbodyJM(hmdat(:,1))
      M    = ft_headcoordinates(hmdat(1:3,1),hmdat(4:6,1),hmdat(7:9,1));
      grad = ft_transform_geometry(M, data.grad);
    end

    volcfg.grad = grad;
    %compute volume conductor
    [volold, grad] = prepare_headmodel(volcfg);
    %compute forward model
    lfold = ft_compute_leadfield(pos, grad, volold);
    %compute projection matrix
    [realign, noalign, bkalign] = computeprojection(lfold, lfnew, cfg.pruneratio, cfg.verify);
  end
  data.realign{i} = realign * data.trial{i};
  if strcmp(cfg.verify, 'yes')
    % also compute the residual variance when interpolating
    [id,it]   = match_str(data.grad.label, template.grad.label);
    rvrealign = rv(data.trial{i}(id,:), data.realign{i}(it,:));
    fprintf('original -> template             RV %.2f %%\n', 100 * mean(rvrealign));
    datnoalign = noalign * data.trial{i};
    datbkalign = bkalign * data.trial{i};
    rvnoalign = rv(data.trial{i}, datnoalign);
    rvbkalign = rv(data.trial{i}, datbkalign);
    fprintf('original             -> original RV %.2f %%\n', 100 * mean(rvnoalign));
    fprintf('original -> template -> original RV %.2f %%\n', 100 * mean(rvbkalign));
  end
end

% plot the topography before and after the realignment
if strcmp(cfg.feedback, 'yes')

  ft_warning('showing MEG topography (RMS value over time) in the first trial only');
  Nchan = length(data.grad.label);
  [id,it]   = match_str(data.grad.label, template.grad.label);
  pos1 = data.grad.chanpos(id,:);
  pos2 = template.grad.chanpos(it,:);
  prj1 = elproj(pos1); tri1 = delaunay(prj1(:,1), prj1(:,2));
  prj2 = elproj(pos2); tri2 = delaunay(prj2(:,1), prj2(:,2));

  switch cfg.topoparam
    case 'rms'
      p1 = sqrt(mean(data.trial{1}(id,:).^2, 2));
      p2 = sqrt(mean(data.realign{1}(it,:).^2, 2));
    case 'svd'
      [u, s, v] = svd(data.trial{1}(id,:)); p1 = u(:,1);
      [u, s, v] = svd(data.realign{1}(it,:)); p2 = u(:,1);
    otherwise
      ft_error('unsupported cfg.topoparam');
  end

  X = [pos1(:,1) pos2(:,1)]';
  Y = [pos1(:,2) pos2(:,2)]';
  Z = [pos1(:,3) pos2(:,3)]';

  % show figure with old an new helmets, volume model and source positions
  figure
  hold on
  ft_plot_headmodel(volold);
  plot3(sourcemodel.pos(:,1),sourcemodel.pos(:,2),sourcemodel.pos(:,3),'b.');
  plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'r.') % original positions
  plot3(pos2(:,1), pos2(:,2), pos2(:,3), 'g.') % template positions
  line(X,Y,Z, 'color', 'black');
  view(-90, 90);

  % show figure with data on old helmet location
  figure
  hold on
  plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'r.') % original positions
  plot3(pos2(:,1), pos2(:,2), pos2(:,3), 'g.') % template positions
  line(X,Y,Z, 'color', 'black');
  axis equal; axis vis3d
  bnd1 = [];
  bnd1.pos = pos1;
  bnd1.tri = tri1;
  ft_plot_mesh(bnd1,'vertexcolor',p1,'edgecolor','none')
  title('RMS, before realignment')
  view(-90, 90)

  % show figure with data on new helmet location
  figure
  hold on
  plot3(pos1(:,1), pos1(:,2), pos1(:,3), 'r.') % original positions
  plot3(pos2(:,1), pos2(:,2), pos2(:,3), 'g.') % template positions
  line(X,Y,Z, 'color', 'black');
  axis equal; axis vis3d
  bnd2 = [];
  bnd2.pos = pos2;
  bnd2.tri = tri2;
  ft_plot_mesh(bnd2,'vertexcolor',p2,'edgecolor','none')
  title('RMS, after realignment')
  view(-90, 90)
end

% store the realigned data in a new structure
interp.label   = template.grad.label;
interp.grad    = template.grad;   % replace with the template gradiometer array
interp.trial   = data.realign;    % remember the processed data
interp.fsample = data.fsample;
interp.time    = data.time;

% add the rest channels back to the data, these were not interpolated
if ~isempty(rest.label)
  fprintf('adding %d non-MEG channels back to the data (', length(rest.label));
  fprintf('%s, ', rest.label{1:end-1});
  fprintf('%s)\n', rest.label{end});
  for trial=1:length(rest.trial)
    interp.trial{trial} = [interp.trial{trial}; rest.trial{trial}];
  end
  interp.label = [interp.label; rest.label];
end

% copy the trial specific information into the output
if isfield(data, 'trialinfo')
  interp.trialinfo = data.trialinfo;
end

% copy the sampleinfo field as well
if isfield(data, 'sampleinfo')
  interp.sampleinfo = data.sampleinfo;
end

% convert back to input type if necessary
switch dtype
  case 'timelock'
    interp = ft_checkdata(interp, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = interp;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction that computes the projection matrix(ces)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [realign, noalign, bkalign] = computeprojection(lfold, lfnew, pruneratio, verify)

% compute this inverse only once, although it is used twice
tmp = prunedinv(lfold, pruneratio);
% compute the three interpolation matrices
fprintf('computing interpolation matrix #1\n');
realign = lfnew * tmp;
if strcmp(verify, 'yes')
  fprintf('computing interpolation matrix #2\n');
  noalign = lfold * tmp;
  fprintf('computing interpolation matrix #3\n');
  bkalign = lfold * prunedinv(lfnew, pruneratio) * realign;
else
  noalign = [];
  bkalign = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction that computes the inverse using a pruned SVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lfi] = prunedinv(lf, r)
[u, s, v] = svd(lf);
if r<1
  % treat r as a ratio
  p = find(s<(s(1,1)*r) & s~=0);
else
  % treat r as the number of spatial components to keep
  diagels = 1:(min(size(s))+1):(min(size(s)).^2);
  p       = diagels((r+1):end);
end
fprintf('pruning %d from %d, i.e. removing the %d smallest spatial components\n', length(p), min(size(s)), length(p));
s(p) = 0;
s(find(s~=0)) = 1./s(find(s~=0));
lfi = v * s' * u';
