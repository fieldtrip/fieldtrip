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
% You should specify the volume conductor model, see FT_FETCH_VOL, and the 
% sensor information, see FT_FETCH_SENS.
%
% If you start with a grid search, you should specify the grid locations at
% which a test dipole will be placed. The positions of the dipoles can be
% specified as a regular 3-D grid that is aligned with the axes of the head
% coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
% Alternatively a complete grid with dipole positions and precomputed
% leadfields can be specified
%   cfg.grid            = structure, see FT_PREPARE_LEADFIELD
% or the position of a few dipoles at locations of interest can be
% specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = Nx3 matrix with position of each source
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%   cfg.grid.inside     = vector with indices of the sources inside the brain (optional)
%   cfg.grid.outside    = vector with indices of the sources outside the brain (optional)
%
% If you do not start with a grid search, you have to give a starting location
% for the nonlinear search
%   cfg.dip.pos     = initial dipole position, matrix of Ndipoles x 3
%
% The conventional approach is to fit dipoles to event-related averages, which
% within fieldtrip can be obtained from the FT_TIMELOCKANALYSIS or from
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
%   cfg.dipfit.display  = level of display, can be 'off', 'iter', 'notify' or 'final' (default = 'iter')
%   cfg.dipfit.optimfun = function to use, can be 'fminsearch' or 'fminunc' (default is determined automatic)
%   cfg.dipfit.maxiter  = maximum number of function evaluations allowed (default depends on the optimfun)
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_LEADFIELD, FT_FETCH_SENS

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
%
% This function depends on FT_PREPARE_DIPOLE_GRID which has the following options:
% cfg.grid.xgrid (default set in FT_PREPARE_DIPOLE_GRID: cfg.grid.xgrid = 'auto'), documented
% cfg.grid.ygrid (default set in FT_PREPARE_DIPOLE_GRID: cfg.grid.ygrid = 'auto'), documented
% cfg.grid.zgrid (default set in FT_PREPARE_DIPOLE_GRID: cfg.grid.zgrid = 'auto'), documented
% cfg.grid.resolution, documented
% cfg.grid.pos, documented
% cfg.grid.dim, documented
% cfg.grid.inside, documented
% cfg.grid.outside, documented
% cfg.symmetry, documented
% cfg.mri
% cfg.smooth
% cfg.sourceunits
% cfg.threshold
%
% This function depends on FT_PREPARE_VOL_SENS which has the following options:
% cfg.channel (default set in FT_DIPOLEFITTING: cfg.channel = 'all'), documented
% cfg.elec, documented
% cfg.elecfile, documented
% cfg.grad, documented
% cfg.gradfile, documented
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented

% Copyright (C) 2004-2006, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar data

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'timelock', 'freq', 'comp'}, 'feedback', 'yes');

% set the defaults
if ~isfield(cfg, 'channel'),     cfg.channel = 'all';        end
if ~isfield(cfg, 'component'),   cfg.component = [];         end  % for comp input
if ~isfield(cfg, 'frequency'),   cfg.frequency = [];         end  % for freq input
if ~isfield(cfg, 'latency'),     cfg.latency = 'all';        end
if ~isfield(cfg, 'feedback'),    cfg.feedback = 'text';      end
if ~isfield(cfg, 'gridsearch'),  cfg.gridsearch = 'yes';     end
if ~isfield(cfg, 'nonlinear'),   cfg.nonlinear = 'yes';      end
if ~isfield(cfg, 'symmetry'),    cfg.symmetry = [];          end

% put the low-level options pertaining to the dipole grid (used for initial scanning) in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'grid'});

% the default for this depends on the data type
if ~isfield(cfg, 'model'),
  if ~isempty(cfg.component)
    % each component is fitted independently
    cfg.model = 'moving';
  elseif ~isempty(cfg.latency)
    % fit the data with a dipole at one location
    cfg.model = 'regional';
  elseif ~isempty(cfg.latency)
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
    error('symmetry constraints are only supported for two-dipole models');
  elseif strcmp(cfg.symmetry, 'x')
    % this structure is passed onto the low-level ft_dipole_fit function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 -1 1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 -x1 y1 z1]
  elseif strcmp(cfg.symmetry, 'y')
    % this structure is passed onto the low-level ft_dipole_fit function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 1 -1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 -y1 z1]
  elseif strcmp(cfg.symmetry, 'z')
    % this structure is passed onto the low-level ft_dipole_fit function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 1 1 -1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 y1 -z1]
  else
    error('unrecognized symmetry constraint');
  end
elseif ~isfield(cfg, 'dipfit') || ~isfield(cfg.dipfit, 'constr')
  % no symmetry constraints have been specified
  cfg.dipfit.constr = [];
end

if isfield(data, 'topolabel')
  % this looks like a component analysis
  iscomp = 1;
  % transform the data into a representation on which the timelocked dipole fit can perform its trick
  data = comp2timelock(cfg, data);
else
  iscomp = 0;
end

if isfield(data, 'freq')
  % this looks like a frequency analysis
  isfreq = 1;
  % transform the data into a representation on which the timelocked dipole fit can perform its trick
  data = freq2timelock(cfg, data);
else
  isfreq = 0;
end

% prepare the volume conduction model and the sensor array
% this updates the configuration with the appropriate fields
[vol, sens, cfg] = prepare_headmodel(cfg, data);

% select the desired channels, the order should be the same as in the sensor structure
[selsens, seldata] = match_str(sens.label, data.label);
Vdata = data.avg(seldata, :);

if iscomp
  % select the desired component topographies
  Vdata = Vdata(:, cfg.component);
elseif isfreq
  % the desired frequency topographies have already been selected
  Vdata = Vdata(:, :);
else
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
fprintf('selected %d channels\n', nchans);
fprintf('selected %d topographies\n', ntime);

if nchans<cfg.numdipoles*3
  warning('not enough channels to perform a dipole fit');
end

if ntime<1
  error('no spatial topography selected');
end

% check whether EEG is average referenced
if ft_senstype(sens, 'eeg')
  if any(rv(Vdata, avgref(Vdata))>0.001)
    warning('the EEG data is not average referenced, correcting this');
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
  error('inconsistent number of dipoles in configuration')
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
  else
    error('dipole scanning is only possible for a single dipole or a symmetric dipole pair');
  end

  % construct the dipole grid on which the gridsearch will be done
  tmpcfg = [];
  tmpcfg.vol  = vol;
  tmpcfg.grad = sens; % this can be electrodes or gradiometers
  % copy all options that are potentially used in ft_prepare_sourcemodel
  try, tmpcfg.grid        = cfg.grid;         end
  try, tmpcfg.mri         = cfg.mri;          end
  try, tmpcfg.headshape   = cfg.headshape;    end
  try, tmpcfg.tightgrid   = cfg.tightgrid;    end
  try, tmpcfg.symmetry    = cfg.symmetry;     end
  try, tmpcfg.smooth      = cfg.smooth;       end
  try, tmpcfg.threshold   = cfg.threshold;    end
  try, tmpcfg.spheremesh  = cfg.spheremesh;   end
  try, tmpcfg.inwardshift = cfg.inwardshift;  end
  try, tmpcfg.sourceunits = cfg.sourceunits;  end
  grid = ft_prepare_sourcemodel(tmpcfg);

  ft_progress('init', cfg.feedback, 'scanning grid');
  for i=1:length(grid.inside)
    ft_progress(i/length(grid.inside), 'scanning grid location %d/%d\n', i, length(grid.inside));
    indx = grid.inside(i);
    if isfield(grid, 'leadfield')
      % reuse the previously computed leadfield
      lf = grid.leadfield{indx};
    else
      lf = ft_compute_leadfield(grid.pos(indx,:), sens, vol);
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
        grid.error(indx,1) = sum(sum(((eye(nchans)-lf*pinv(lf))*Vdata).^2));
      case 'moving'
        % remember the error for each latency independently
        grid.error(indx,:) = sum(((eye(nchans)-lf*pinv(lf))*Vdata).^2);
      otherwise
        error('unsupported cfg.model');
    end % switch model
  end % looping over the grid
  ft_progress('close');
  
  switch cfg.model
    case 'regional'
      % find the grid point(s) with the minimum error
      [err, indx] = min(grid.error(grid.inside));
      dip.pos = grid.pos(grid.inside(indx),:);          % note that for a symmetric dipole pair this results in a vector
      dip.pos = reshape(dip.pos,3,cfg.numdipoles)';    % convert to a Nx3 array
      dip.mom = zeros(cfg.numdipoles*3,1);              % set the dipole moment to zero
      if cfg.numdipoles==1
        fprintf('found minimum after scanning on grid point [%g %g %g]\n', dip.pos(1), dip.pos(2), dip.pos(3));
      elseif cfg.numdipoles==2
        fprintf('found minimum after scanning on grid point [%g %g %g; %g %g %g]\n', dip.pos(1), dip.pos(2), dip.pos(3), dip.pos(4), dip.pos(5), dip.pos(6));
      end
    case 'moving'
      for t=1:ntime
        % find the grid point(s) with the minimum error
        [err, indx] = min(grid.error(grid.inside,t));
        dip(t).pos = grid.pos(grid.inside(indx),:);           % note that for a symmetric dipole pair this results in a vector
        dip(t).pos = reshape(dip(t).pos,3,cfg.numdipoles)';  % convert to a Nx3 array
        dip(t).mom = zeros(cfg.numdipoles*3,1);               % set the dipole moment to zero
        if cfg.numdipoles==1
          fprintf('found minimum after scanning for topography %d on grid point [%g %g %g]\n', t, dip(t).pos(1), dip(t).pos(2), dip(t).pos(3));
        elseif cfg.numdipoles==2
          fprintf('found minimum after scanning for topography %d on grid point [%g %g %g; %g %g %g]\n', t, dip(t).pos(1), dip(t).pos(2), dip(t).pos(3), dip(t).pos(4), dip(t).pos(5), dip(t).pos(6));
        end
      end
    otherwise
      error('unsupported cfg.model');
  end % switch model
end % if gridsearch

if strcmp(cfg.gridsearch, 'no')
  % use the initial guess supplied in the configuration for the remainder
  switch cfg.model
    case 'regional'
      dip = struct(cfg.dip);
    case 'moving'
      for t=1:ntime
        dip(t) = struct(cfg.dip);
      end
    otherwise
      error('unsupported cfg.model');
  end % switch model
end

% multiple dipoles can be represented either as a 1x(N*3) vector or as a  Nx3 matrix,
% i.e. [x1 y1 z1 x2 y2 z2] or [x1 y1 z1; x2 y2 z2]
switch cfg.model
  case 'regional'
    dip = fixdipole(dip);
  case 'moving'
    for t=1:ntime
      dip(t) = fixdipole(dip(t));
    end
  otherwise
    error('unsupported cfg.model');
end % switch model

if isfield(cfg, 'dipfit')
  % convert the structure with the additional low-level options into key-value pairs
  optarg = ft_cfg2keyval(getfield(cfg, 'dipfit'));
else
  % no additional low-level options were specified
  optarg = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the non-linear fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(cfg.nonlinear, 'yes')
  switch cfg.model
    case 'regional'
      % perform the non-linear dipole fit for all latencies together
      % catch errors due to non-convergence
      try
        dip = dipole_fit(dip, sens, vol, Vdata, optarg{:});
        success = 1;
        if cfg.numdipoles==1
          fprintf('found minimum after non-linear optimization on [%g %g %g]\n', dip.pos(1), dip.pos(2), dip.pos(3));
        elseif cfg.numdipoles==2
          fprintf('found minimum after non-linear optimization on [%g %g %g; %g %g %g]\n', dip.pos(1,1), dip.pos(1,2), dip.pos(1,3), dip.pos(2,1), dip.pos(2,2), dip.pos(2,3));
        end
      catch
        success = 0;
        disp(lasterr);
      end
    case 'moving'
      % perform the non-linear dipole fit for each latency independently
      % instead of using dip(t) =  ft_dipole_fit(dip(t),...), I am using temporary variables dipin and dipout
      % this is to prevent errors of the type "Subscripted assignment between dissimilar structures"
      dipin = dip;
      for t=1:ntime
        % catch errors due to non-convergence
        try
          dipout(t) = dipole_fit(dipin(t), sens, vol, Vdata(:,t), optarg{:});
          success(t) = 1;
          if cfg.numdipoles==1
            fprintf('found minimum after non-linear optimization for topography %d on [%g %g %g]\n', t, dipout(t).pos(1), dipout(t).pos(2), dipout(t).pos(3));
          elseif cfg.numdipoles==2
            fprintf('found minimum after non-linear optimization for topography %d on [%g %g %g; %g %g %g]\n', t, dipout(t).pos(1,1), dipout(t).pos(1,2), dipout(t).pos(1,3), dipout(t).pos(2,1), dipout(t).pos(2,2), dipout(t).pos(2,3));
          end
        catch
          dipout(t).pos = dipin(t).pos;
          dipout(t).mom = dipin(t).mom;
          success(t) = 0;
          disp(lasterr);
        end
      end
      dip = dipout;
      clear dipin dipout
    otherwise
      error('unsupported cfg.model');
  end % switch model
end % if nonlinear

if strcmp(cfg.nonlinear, 'no')
  % the optimal dipole positions are either obrained from scanning or from the initial configured seed
  switch cfg.model
    case 'regional'
      success = 1;
    case 'moving'
      success = ones(1,ntime);
    otherwise
      error('unsupported cfg.model');
      
  end % switch model
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the model potential distribution and the residual variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.model
  case 'regional'
    if success
      % re-compute the leadfield in order to compute the model potential and dipole moment
      lf = ft_compute_leadfield(dip.pos, sens, vol);
      % compute all details of the final dipole model
      dip.mom = pinv(lf)*Vdata;
      dip.pot = lf*dip.mom;
      dip.rv  = rv(Vdata, dip.pot);
      Vmodel  = dip.pot;
    end
  case 'moving'
    for t=1:ntime
      if success(t)
        % re-compute the leadfield in order to compute the model potential and dipole moment
        lf = ft_compute_leadfield(dip(t).pos, sens, vol);
        % compute all details of the final dipole model
        dip(t).mom = pinv(lf)*Vdata(:,t);
        dip(t).pot = lf*dip(t).mom;
        dip(t).rv  = rv(Vdata(:,t), dip(t).pot);
        Vmodel(:,t) = dip(t).pot;
      end
    end
  otherwise
    error('unsupported cfg.model');
end % switch model

switch cfg.model
  case 'regional'
    if isfreq
      % the matrix with the dipole moment is encrypted and cannot be interpreted straight away
      % reconstruct the frequency representation of the data at the source level
      [dip.pow, dip.csd, dip.fourier] = timelock2freq(dip.mom);
    end
  case 'moving'
    if isfreq
      % although this is technically possible so far, it does not make any sense
      warning('a moving dipole model in the frequency domain is not supported');
    end
  otherwise
    error('unsupported cfg.model');
end % switch model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source.label  = cfg.channel; % these channels were used in fitting
source.dip    = dip;
source.Vdata  = Vdata;  % FIXME this should be renamed (if possible w.r.t. EEGLAB)
source.Vmodel = Vmodel; % FIXME this should be renamed (if possible w.r.t. EEGLAB)

% assign a latency, frequeny or component axis to the output
if iscomp
  source.component = cfg.component;
  % FIXME assign Vdata to an output variable, idem for the model potential
elseif isfreq
  source.freq   = cfg.frequency;
  source.dimord = 'chan_freq';
  % FIXME assign Vdata to an output variable, idem for the model potential
else
  tbeg = nearest(data.time, cfg.latency(1));
  tend = nearest(data.time, cfg.latency(end));
  source.time   = data.time(tbeg:tend);
  source.dimord = 'chan_time';
end

% FIXME why would this be done?
if isfield(data, 'grad')
  % copy the gradiometer array along
  source.grad = data.grad;
end
if isfield(data, 'elec')
  % copy the electrode array along
  source.elec = data.elec;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data
ft_postamble history source
ft_postamble savevar source
