function [source] = ft_dipolefitting(cfg, data)

% FT_DIPOLEFITTING perform grid search and non-linear fit with one or multiple
% dipoles and try to find the location where the dipole model is best able
% to explain the measured EEG or MEG topography.
%
% This function will initially scan the whole brain with a single dipole on
% a regular coarse grid, and subsequently start at the most optimal location
% with a non-linear search. Alternatively you can specify the initial
% location of the dipole(s) and the non-linear search will start from there.
%
% Use as
%   [source] = ft_dipolefitting(cfg, data)
%
% The configuration has the following general fields
%   cfg.numdipoles  = number, default is 1
%   cfg.symmetry    = 'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
%   cfg.gridsearch  = 'yes' or 'no', perform global search for initial
%                     guess for the dipole parameters (default = 'yes')
%   cfg.nonlinear   = 'yes' or 'no', perform nonlinear search for optimal
%                     dipole parameters (default = 'yes')
%
% If you start with a grid search, the complete grid with dipole positions is
% constructed using FT_PREPARE_SOURCEMODEL. It can be specified as as a regular 3-D
% grid that is aligned with the axes of the head coordinate system using
%   cfg.xgrid               = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.ygrid               = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.zgrid               = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.resolution          = number (e.g. 1 cm) for automatic grid generation
% If the source model destribes a triangulated cortical sheet, it is described as
%   cfg.sourcemodel.pos     = N*3 matrix with the vertex positions of the cortical sheet
%   cfg.sourcemodel.tri     = M*3 matrix that describes the triangles connecting the vertices
% Alternatively the position of a few dipoles at locations of interest can be
% user-specified, for example obtained from an anatomical or functional MRI
%   cfg.sourcemodel.pos     = N*3 matrix with position of each source
%   cfg.sourcemodel.inside  = N*1 vector with boolean value whether grid point is inside brain (optional)
%   cfg.sourcemodel.dim     = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%
% If you do not start with a grid search, you have to give a starting location
% for the nonlinear search
%   cfg.dip.pos     = initial dipole position, matrix of Ndipoles x 3
%
% The conventional approach is to fit dipoles to event-related averages, which
% within FieldTrip can be obtained from the FT_TIMELOCKANALYSIS or from
% the FT_TIMELOCKGRANDAVERAGE function. This has the additional options
%   cfg.latency     = [begin end] in seconds or 'all' (default = 'all')
%   cfg.model       = 'moving' or 'regional'
% A moving dipole model has a different position (and orientation) for each
% timepoint, or for each component. A regional dipole model has the same
% position for each timepoint or component, and a different orientation.
%
% You can also fit dipoles to the spatial topographies of an independent
% component analysis, obtained from the FT_COMPONENTANALYSIS function.
% This has the additional options
%   cfg.component   = array with numbers (can be empty -> all)
%
% You can also fit dipoles to the spatial topographies that are present
% in the data in the frequency domain, which can be obtained using the
% FT_FREQANALYSIS function. This has the additional options
%   cfg.frequency   = single number (in Hz)
%
% Low level details of the fitting can be specified in the cfg.dipfit structure
%   cfg.dipfit.display      = level of display, can be 'off', 'iter', 'notify' or 'final' (default = 'iter')
%   cfg.dipfit.optimfun     = function to use, can be 'fminsearch' or 'fminunc' (default is determined automatic)
%   cfg.dipfit.maxiter      = maximum number of function evaluations allowed (default depends on the optimfun)
%   cfg.dipfit.checkinside  = boolean, check that the dipole remains in the source compartment (default = false)
%
% Optionally, you can modify the leadfields by reducing the rank, i.e. remove the weakest orientation
%   cfg.reducerank    = 'no', or number (default = 3 for EEG, 2 for MEG)
%   cfg.backproject   = 'yes' or 'no',  determines when reducerank is applied whether the
%                       lower rank leadfield is projected back onto the original linear
%                       subspace, or not (default = 'yes')
%
% The volume conduction model of the head should be specified as
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition or filename, see FT_READ_SENS
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_LEADFIELD, FT_PREPARE_HEADMODEL

% TODO change the output format, more suitable would be something like:
% dip.label
% dip.time
% dip.avg (instead of Vdata)
% dip.dip.pos
% dip.dip.mom
% dip.dip.model, or dip.dip.avg
% dip.dimord

% Undocumented local options:
%   cfg.dipfit.constr   = Source model constraints, depends on cfg.symmetry
% Optionally, you can include a noise covariance structure to sphere the data (is useful when using both
% magnetometers and gradiometers to fit your dipole)
%   cfg.dipfit.noisecov       = noise covariance matrix, see e.g. FT_TIMELOCK_ANALYSIS

% Copyright (C) 2004-2013, Robert Oostenveld
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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'comp', 'timelock', 'freq'}, 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamed',    {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed',    {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed',    {'optofile', 'opto'});
cfg = ft_checkconfig(cfg, 'renamed',    {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed',    {'grid',    'sourcemodel'});

% get the defaults
cfg.channel         = ft_getopt(cfg, 'channel', 'all');
cfg.component       = ft_getopt(cfg, 'component', 'all');   % for comp input
cfg.frequency       = ft_getopt(cfg, 'frequency');          % for freq input
cfg.latency         = ft_getopt(cfg, 'latency', 'all');     % for timelock input
cfg.feedback        = ft_getopt(cfg, 'feedback', 'text');
cfg.gridsearch      = ft_getopt(cfg, 'gridsearch', 'yes');
cfg.nonlinear       = ft_getopt(cfg, 'nonlinear', 'yes');
cfg.symmetry        = ft_getopt(cfg, 'symmetry');
cfg.dipfit          = ft_getopt(cfg, 'dipfit', []);     % the default for this is handled below

cfg = ft_checkconfig(cfg, 'renamed',    {'tightgrid', 'tight'});  % this is moved to cfg.sourcemodel.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed',    {'sourceunits', 'unit'}); % this is moved to cfg.sourcemodel.unit  by the subsequent createsubcfg

% put the low-level options pertaining to the sourcemodel in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg', {'sourcemodel'});
% move some fields from cfg.sourcemodel back to the top-level configuration
cfg = ft_checkconfig(cfg, 'createtopcfg', {'sourcemodel'});

% determine data type
iscomp = ft_datatype(data, 'comp');           % it can also be raw+comp, timelock+comp or freq+comp
isfreq = ft_datatype(data, 'freq');           % it might also be freq+comp, in that case it should be treated as component data
istimelock = ft_datatype(data, 'timelock');   % it might also be timelock+comp, in that case it should be treated as component data

% the default for this depends on the data type
if ~isfield(cfg, 'model')
  if iscomp
    % each component is fitted independently
    cfg.model = 'moving';
  elseif isfreq
    % fit the data with a dipole at one location
    cfg.model = 'regional';
  elseif istimelock
    % fit the data with a dipole at one location
    cfg.model = 'regional';
  end
end

if ~isfield(cfg, 'numdipoles')
  if isfield(cfg, 'dip')
    cfg.numdipoles = size(cfg.dip(1).pos,1);
  else
    cfg.numdipoles = 1;
  end
end

% set up the symmetry constraints
if ~isempty(cfg.symmetry)
  if cfg.numdipoles~=2
    ft_error('symmetry constraints are only supported for two-dipole models');
  elseif strcmp(cfg.symmetry, 'x')
    % this structure is passed onto the low-level FT_INVERSE_DIPOLEFIT function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 -1 1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 -x1 y1 z1]
  elseif strcmp(cfg.symmetry, 'y')
    % this structure is passed onto the low-level FT_INVERSE_DIPOLEFIT function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 1 -1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 -y1 z1]
  elseif strcmp(cfg.symmetry, 'z')
    % this structure is passed onto the low-level FT_INVERSE_DIPOLEFIT function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 1 1 -1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 y1 -z1]
  else
    ft_error('unrecognized symmetry constraint');
  end
elseif ~isfield(cfg, 'dipfit') || ~isfield(cfg.dipfit, 'constr')
  % no symmetry constraints have been specified
  cfg.dipfit.constr = [];
end

if ft_getopt(cfg.dipfit.constr, 'sequential', false) && strcmp(cfg.model, 'moving')
  ft_error('the moving dipole model does not combine with the sequential constraint')
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3119
end

if iscomp
  % transform the data into a representation on which the timelocked dipole fit can perform its trick
  data = comp2timelock(cfg, data);
  
  % default component selection is all components
  if ischar(cfg.component) && strcmp(cfg.component, 'all')
    cfg.component = (1:size(data.avg, 2));
  end
elseif isfreq
  % transform the data into a representation on which the timelocked dipole fit can perform its trick
  data = freq2timelock(cfg, data);
elseif istimelock
  % no transformation is needed
end
  
% collect and preprocess the electrodes/gradiometer and head model
% this will also update cfg.channel to match the electrodes/gradiometers
[headmodel, sens, cfg] = prepare_headmodel(cfg, data);

% construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD and FT_INVERSE_DIPOLEFIT
leadfieldopt = {};
leadfieldopt = ft_setopt(leadfieldopt, 'reducerank',     ft_getopt(cfg, 'reducerank'));
leadfieldopt = ft_setopt(leadfieldopt, 'backproject',    ft_getopt(cfg, 'backproject'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(cfg, 'normalize'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(cfg, 'normalizeparam'));
leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(cfg, 'weight'));

% construct the low-level options for the dipole fitting as key-value pairs, these are passed to FT_INVERSE_DIPOLEFIT
dipfitopt = ft_cfg2keyval(cfg.dipfit);

% select the desired channels, ordered according to the sensor structure or configuration
[selcfg, seldata] = match_str(cfg.channel, data.label);
% take the selected channels from the data structure
Vdata = data.avg(seldata, :);

% sphere the date using the noise covariance matrix supplied, if any
% this affects both the gridsearch and the nonlinear optimization
noisecov = ft_getopt(cfg.dipfit, 'noisecov');
if ~isempty(noisecov)
  [u, s] = svd(noisecov);
  tol = max(size(noisecov)) * eps(norm(s, inf));
  s = diag(s);
  r1 = sum(s > tol) + 1;
  s(1:(r1 - 1)) = 1 ./ sqrt(s(1:(r1 - 1)));
  s(r1:end)     = 0;
  sphere = diag(s) * u';
  % apply the sphering to the data
  Vdata = sphere * Vdata;
  % apply the sphering as a pre-multiplication to the sensor definition
  montage = [];
  montage.labelold = cfg.channel;
  montage.labelnew = cfg.channel;
  montage.tra = sphere;
  sens = ft_apply_montage(sens, montage, 'balancename', 'sphering');
end

if iscomp
  % select the desired component topographies
  Vdata = Vdata(:, cfg.component);
elseif isfreq
  % the desired frequencies have already been selected
  Vdata = Vdata(:, :);
elseif istimelock
  % select the desired latencies
  if ischar(cfg.latency) && strcmp(cfg.latency, 'all')
    cfg.latency = data.time([1 end]);
  end
  tbeg = nearest(data.time, cfg.latency(1));
  tend = nearest(data.time, cfg.latency(end));
  cfg.latency = [data.time(tbeg) data.time(tend)];
  Vdata = Vdata(:, tbeg:tend);
end

nchans = size(Vdata,1);
ntime  = size(Vdata,2);
Vmodel = zeros(nchans, ntime);
ft_info('selected %d channels\n', nchans);
ft_info('selected %d topographies\n', ntime);

if nchans<cfg.numdipoles*3
  ft_warning('not enough channels to perform a dipole fit');
end

if ntime<1
  ft_error('no spatial topography selected');
end

% check whether EEG is average referenced
if ft_senstype(sens, 'eeg')
  if any(rv(Vdata, avgref(Vdata))>0.001)
    ft_warning('the EEG data is not average referenced, correcting this');
  end
  Vdata = avgref(Vdata);
end

% set to zeros if no initial dipole was specified
if ~isfield(cfg, 'dip')
  cfg.dip.pos = zeros(cfg.numdipoles, 3);
  cfg.dip.mom = zeros(3*cfg.numdipoles, 1);
end

% set to zeros if no initial dipole position was specified
if ~isfield(cfg.dip, 'pos')
  cfg.dip.pos = zeros(cfg.numdipoles, 3);
end

% set to zeros if no initial dipole moment was specified
if ~isfield(cfg.dip, 'mom')
  cfg.dip.mom = zeros(3*cfg.numdipoles, 1);
end

% check the specified dipole model
if numel(cfg.dip.pos)~=cfg.numdipoles*3 || numel(cfg.dip.mom)~=cfg.numdipoles*3
  ft_error('inconsistent number of dipoles in configuration')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the dipole scan, this is usefull for generating an initial guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(cfg.gridsearch, 'yes')
  % test whether we have a valid configuration for dipole scanning
  if cfg.numdipoles==1
    % this is ok
  elseif cfg.numdipoles==2 && ~isempty(cfg.dipfit.constr)
    % this is also ok
  elseif isfield(cfg.sourcemodel, 'pos') && size(cfg.sourcemodel.pos,2)==cfg.numdipoles*3
    % this is also ok
  else
    ft_error('dipole scanning is only possible for a single dipole or a symmetric dipole pair');
  end
  
  if isfield(cfg.sourcemodel, 'leadfield')
    ft_notice('using precomputed leadfields for the gridsearch');

    sourcemodel = keepfields(cfg.sourcemodel, {'pos', 'tri', 'dim', 'inside', 'leadfield', 'leadfielddimord', 'label'});
    
    % select the channels corresponding to the data and the user configuration
    tmpcfg = keepfields(cfg, 'channel');
    sourcemodel = ft_selectdata(tmpcfg, sourcemodel);
    
    % sort the channels to be consistent with the data
    [dum, chansel] = match_str(data.label, sourcemodel.label);
    sourcemodel.label = sourcemodel.label(chansel);
    for i=1:numel(sourcemodel.leadfield)
      if ~isempty(sourcemodel.leadfield{i})
        sourcemodel.leadfield{i} = sourcemodel.leadfield{i}(chansel, :);
      end
    end
    
    % ensure that the channels are consistent with the data
    assert(isequal(sourcemodel.label, cfg.channel), 'cannot match the channels in the sourcemodel to those in the data')
    
  else
    ft_notice('computing the leadfields for the gridsearch on the fly');
    
    % construct the dipole positions on which the source reconstruction will be done
    tmpcfg           = keepfields(cfg, {'sourcemodel', 'mri', 'headshape', 'symmetry', 'smooth', 'threshold', 'spheremesh', 'inwardshift', 'xgrid' 'ygrid', 'zgrid', 'resolution', 'tight', 'warpmni', 'template', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
    tmpcfg.headmodel = headmodel;
    if ft_senstype(sens, 'eeg')
      tmpcfg.elec = sens;
    elseif ft_senstype(sens, 'meg')
      tmpcfg.grad = sens;
    end
    sourcemodel = ft_prepare_sourcemodel(tmpcfg);
    
  end % if precomputed leadfield or not

  ngrid = size(sourcemodel.pos,1);
  
  switch cfg.model
    case 'regional'
      sourcemodel.error = nan(ngrid, 1);
    case 'moving'
      sourcemodel.error = nan(ngrid, ntime);
    otherwise
      ft_error('unsupported cfg.model');
  end
  
  insideindx = find(sourcemodel.inside);
  ft_progress('init', cfg.feedback, 'scanning grid');
  for i=1:length(insideindx)
    ft_progress(i/length(insideindx), 'scanning grid location %d/%d\n', i, length(insideindx));
    thisindx = insideindx(i);
    if isfield(sourcemodel, 'leadfield')
      % reuse the previously computed leadfield
      lf = sourcemodel.leadfield{thisindx};
    else
      lf = ft_compute_leadfield(sourcemodel.pos(thisindx,:), sens, headmodel, leadfieldopt{:});
    end
    % the model is V=lf*mom+noise, therefore mom=pinv(lf)*V estimates the
    % dipole moment this makes the model potential U=lf*pinv(lf)*V and the
    % model error is norm(V-U) = norm(V-lf*pinv(lf)*V) = norm((eye-lf*pinv(lf))*V)
    if any(isnan(lf(:)))
      % this might happen if one of the dipole locations of the grid is
      % outside the brain compartment
      lf(:) = 0;
    end
    switch cfg.model
      case 'regional'
        % sum the error over all latencies
        sourcemodel.error(thisindx,1) = sum(sum(((eye(nchans)-lf*pinv(lf))*Vdata).^2));
      case 'moving'
        % remember the error for each latency independently
        sourcemodel.error(thisindx,:) = sum(((eye(nchans)-lf*pinv(lf))*Vdata).^2);
      otherwise
        ft_error('unsupported cfg.model');
    end % switch model
  end % looping over the grid
  ft_progress('close');
  
  switch cfg.model
    case 'regional'
      % find the source position with the minimum error
      [err, indx] = min(sourcemodel.error);
      dip.pos = sourcemodel.pos(indx,:);                % note that for a symmetric dipole pair this results in a vector
      dip.pos = reshape(dip.pos,3,cfg.numdipoles)';     % convert to a Nx3 array
      dip.mom = zeros(cfg.numdipoles*3,1);              % set the dipole moment to zero
      if cfg.numdipoles==1
        ft_info('found minimum after scanning on grid point [%g %g %g]\n', dip.pos(1), dip.pos(2), dip.pos(3));
      elseif cfg.numdipoles==2
        ft_info('found minimum after scanning on grid point [%g %g %g; %g %g %g]\n', dip.pos(1,1), dip.pos(1,2), dip.pos(1,3), dip.pos(2,1), dip.pos(2,2), dip.pos(2,3));
      end
      
    case 'moving'
      for t=1:ntime
        % find the source position with the minimum error
        [err, indx] = min(sourcemodel.error(:,t));
        dip(t).pos = sourcemodel.pos(indx,:);                 % note that for a symmetric dipole pair this results in a vector
        dip(t).pos = reshape(dip(t).pos,3,cfg.numdipoles)';   % convert to a Nx3 array
        dip(t).mom = zeros(cfg.numdipoles*3,1);               % set the dipole moment to zero
        if cfg.numdipoles==1
          ft_info('found minimum after scanning for topography %d on grid point [%g %g %g]\n', t, dip(t).pos(1), dip(t).pos(2), dip(t).pos(3));
        elseif cfg.numdipoles==2
          ft_info('found minimum after scanning for topography %d on grid point [%g %g %g; %g %g %g]\n', t, dip(t).pos(1,1), dip(t).pos(1,2), dip(t).pos(1,3), dip(t).pos(2,1), dip(t).pos(2,2), dip(t).pos(2,3));
        end
      end
      
    otherwise
      ft_error('unsupported cfg.model');
  end % switch model
  
elseif strcmp(cfg.gridsearch, 'no')
  % there is no grid needed for dipole scanning
  sourcemodel = [];
  % use the initial guess supplied in the configuration for the remainder
  switch cfg.model
    case 'regional'
      dip = cfg.dip;
    case 'moving'
      for t=1:ntime
        dip(t) = cfg.dip;
      end
    otherwise
      ft_error('unsupported cfg.model');
  end % switch model
  
end % if gridsearch yes/no

% multiple dipoles can be represented either as a 1x(N*3) vector or as a Nx3 matrix,
% i.e. [x1 y1 z1 x2 y2 z2] or [x1 y1 z1; x2 y2 z2]
switch cfg.model
  case 'regional'
    dip = fixdipole(dip);
  case 'moving'
    for t=1:ntime
      dip(t) = fixdipole(dip(t));
    end
  otherwise
    ft_error('unsupported cfg.model');
end % switch model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the non-linear fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(cfg.nonlinear, 'yes')
  switch cfg.model
    case 'regional'
      % perform the non-linear dipole fit for all latencies together
      % catch errors due to non-convergence
      try
        dip = ft_inverse_dipolefit(dip, sens, headmodel, Vdata, dipfitopt{:}, leadfieldopt{:});
        success = 1;
        if cfg.numdipoles==1
          ft_info('found minimum after non-linear optimization on [%g %g %g]\n', dip.pos(1), dip.pos(2), dip.pos(3));
        elseif cfg.numdipoles==2
          ft_info('found minimum after non-linear optimization on [%g %g %g; %g %g %g]\n', dip.pos(1,1), dip.pos(1,2), dip.pos(1,3), dip.pos(2,1), dip.pos(2,2), dip.pos(2,3));
        end
      catch
        success = 0;
        disp(lasterr);
      end
      
    case 'moving'
      % perform the non-linear dipole fit for each latency independently
      % instead of using dip(t) = ft_inverse_dipolefit(dip(t),...), I am using temporary variables dipin and dipout
      % to prevent errors like "Subscripted assignment between dissimilar structures"
      dipin = dip;
      for t=1:ntime
        % catch errors due to non-convergence
        try
          dipout(t) = ft_inverse_dipolefit(dipin(t), sens, headmodel, Vdata(:,t), dipfitopt{:}, leadfieldopt{:});
          success(t) = 1;
          if cfg.numdipoles==1
            ft_info('found minimum after non-linear optimization for topography %d on [%g %g %g]\n', t, dipout(t).pos(1), dipout(t).pos(2), dipout(t).pos(3));
          elseif cfg.numdipoles==2
            ft_info('found minimum after non-linear optimization for topography %d on [%g %g %g; %g %g %g]\n', t, dipout(t).pos(1,1), dipout(t).pos(1,2), dipout(t).pos(1,3), dipout(t).pos(2,1), dipout(t).pos(2,2), dipout(t).pos(2,3));
          end
        catch
          % keep the position and moment according to the initial guess
          dipout(t).pos = dipin(t).pos;
          dipout(t).mom = dipin(t).mom;
          success(t) = 0;
          disp(lasterr);
        end
      end
      dip = dipout;
      clear dipin dipout
    otherwise
      ft_error('unsupported cfg.model');
  end % switch model
end % if nonlinear

if strcmp(cfg.nonlinear, 'no')
  % the optimal dipole positions are either obtained from scanning
  % or from the initial configured specified by the user
  switch cfg.model
    case 'regional'
      success = 1;
    case 'moving'
      success = ones(1,ntime);
    otherwise
      ft_error('unsupported cfg.model');
      
  end % switch model
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the model potential distribution and the residual variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.model
  case 'regional'
    if success
      % re-compute the leadfield in order to compute the model potential and dipole moment
      lf = ft_compute_leadfield(dip.pos, sens, headmodel, leadfieldopt{:});
      if isfield(dip, 'mom') && isfield(dip, 'ampl')
        % the orientation and amplitude have already been estimated, this applies to the case of a fixed dipole orientation
        dip.pot = (lf * dip.mom) * dip.ampl;
      else
        % compute all details of the final dipole model using linear estimation
        dip.mom = pinv(lf)*Vdata;
        dip.pot = lf*dip.mom;
      end
      dip.rv  = rv(Vdata, dip.pot);
      Vmodel  = dip.pot;
    end
  case 'moving'
    for t=1:ntime
      if success(t)
        % re-compute the leadfield in order to compute the model potential and dipole moment
        lf = ft_compute_leadfield(dip(t).pos, sens, headmodel, leadfieldopt{:});
        % compute all details of the final dipole model
        dip(t).mom = pinv(lf)*Vdata(:,t);
        dip(t).pot = lf*dip(t).mom;
        dip(t).rv  = rv(Vdata(:,t), dip(t).pot);
        Vmodel(:,t) = dip(t).pot;
      end
    end
  otherwise
    ft_error('unsupported cfg.model');
end % switch model

switch cfg.model
  case 'regional'
    if isfreq
      % the matrix with the dipole moment is encrypted and cannot be interpreted straight away
      % reconstruct the frequency representation of the data at the source level
      if isfield(dip, 'mom') && isfield(dip, 'ampl')
        % this applies to the case of a fixed dipole orientation
        [dip.pow, dip.csd, dip.fourier] = timelock2freq(dip.mom * dip.ampl);
      else
        [dip.pow, dip.csd, dip.fourier] = timelock2freq(dip.mom);
      end
    end
  case 'moving'
    if isfreq
      % although this is technically possible so far, it does not make any sense
      ft_warning('a moving dipole model in the frequency domain is not supported');
    end
  otherwise
    ft_error('unsupported cfg.model');
end % switch model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source.label  = cfg.channel; % these channels were used in fitting
source.dip    = dip;
source.Vdata  = Vdata;  % FIXME this should be renamed (if possible w.r.t. EEGLAB)
source.Vmodel = Vmodel; % FIXME this should be renamed (if possible w.r.t. EEGLAB)

% the units of the fitted source are the same as the units of the headmodel and the sensor array
for i=1:length(source.dip)
  if isfield(headmodel, 'unit')
    source.dip(i).unit = headmodel.unit;
  elseif isfield(sourcemodel, 'unit')
    source.dip(i).unit = sourcemodel.unit;
  end
end

% assign a latency, frequeny or component axis to the output
if iscomp
  source.component = cfg.component;
  % FIXME assign Vdata to an output variable, idem for the model potential
elseif isfreq
  source.freq   = cfg.frequency;
  source.dimord = 'chan_freq';
  % FIXME assign Vdata to an output variable, idem for the model potential
elseif istimelock
  tbeg = nearest(data.time, cfg.latency(1));
  tend = nearest(data.time, cfg.latency(end));
  source.time   = data.time(tbeg:tend);
  source.dimord = 'chan_time';
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance source
ft_postamble history    source
ft_postamble savevar    source
