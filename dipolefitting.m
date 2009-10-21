function [source] = dipolefitting(cfg, data)

% DIPOLEFITTING perform grid search and non-linear fit with one or multiple
% dipoles and try to find the location where the dipole model is best able
% to explain the measured EEG or MEG topography.
%
% This function will initially scan the whole brain with a single dipole on
% a regular coarse grid, and subsequently start at the most optimal location
% with a non-linear search. Alternatively you can specify the initial
% location of the dipole(s) and the non-linear search will start from there.
%
% Use as
%   [source] = dipolefitting(cfg, data)
%
% The configuration has the following general fields
%   cfg.numdipoles  = number, default is 1
%   cfg.symmetry    = 'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see CHANNELSELECTION for details
%   cfg.gridsearch  = 'yes' or 'no', perform global search for initial
%                     guess for the dipole parameters (default = 'yes')
%   cfg.nonlinear   = 'yes' or 'no', perform nonlinear search for optimal
%                     dipole parameters (default = 'yes')
%
% You should specify the volume conductor model with
%   cfg.hdmfile     = string, file containing the volume conduction model
% or alternatively
%   cfg.vol         = structure with volume conduction model
% If the sensor information is not contained in the data itself you should
% also specify the sensor information using
%   cfg.gradfile    = string, file containing the gradiometer definition
%   cfg.elecfile    = string, file containing the electrode definition
% or alternatively
%   cfg.grad        = structure with gradiometer definition
%   cfg.elec        = structure with electrode definition
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
%   cfg.grid            = structure, see PREPARE_LEADFIELD
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
% within fieldtrip can be obtained from the TIMELOCKANALYSIS or from
% the TIMELOCKGRANDAVERAGE function. This has the additional options
%   cfg.latency     = [begin end] in seconds or 'all' (default = 'all')
%   cfg.model       = 'moving' or 'regional'
% A moving dipole model has a different position (and orientation) for each
% timepoint, or for each component. A regional dipole model has the same
% position for each timepoint or component, and a different orientation.
%
% You can also fit dipoles to the spatial topographies of an independent
% component analysis, obtained from the COMPONENTANALYSIS function.
% This has the additional options
%   cfg.component   = array with numbers (can be empty -> all)
%
% You can also fit dipoles to the spatial topographies that are present
% in the data in the frequency domain, which can be obtained using the
% FREQANALYSIS function. This has the additional options
%   cfg.frequency   = single number (in Hz)
%
% Low level details of the fitting can be specified in the cfg.dipfit structure
%   cfg.dipfit.display  = level of display, can be 'off', 'iter', 'notify' or 'final' (default = 'iter')
%   cfg.dipfit.optimfun = function to use, can be 'fminsearch' or 'fminunc' (default is determined automatic)
%   cfg.dipfit.maxiter  = maximum number of function evaluations allowed (default depends on the optimfun)
%
% See also SOURCEANALYSIS, PREPARE_LEADFIELD

% TODO change the output format, more suitable would be something like:
% dip.label
% dip.time
% dip.avg (instead of Vdata)
% dip.dip.pos
% dip.dip.mom
% dip.dip.model, or dip.dip.avg
% dip.dimord

% Undocumented local options:
% cfg.dipfit.constr   = Source model constraints, depends on cfg.symmetry
%
% This function depends on PREPARE_DIPOLE_GRID which has the following options:
% cfg.grid.xgrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.xgrid = 'auto'), documented
% cfg.grid.ygrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.ygrid = 'auto'), documented
% cfg.grid.zgrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.zgrid = 'auto'), documented
% cfg.grid.resolution, documented
% cfg.grid.pos, documented
% cfg.grid.dim, documented
% cfg.grid.inside, documented
% cfg.grid.outside, documented
% cfg.symmetry, documented
% cfg.mri
% cfg.mriunits
% cfg.smooth
% cfg.sourceunits
% cfg.threshold
%
% This function depends on PREPARE_VOL_SENS which has the following options:
% cfg.channel (default set in DIPOLEFITTING: cfg.channel = 'all'), documented
% cfg.elec, documented
% cfg.elecfile, documented
% cfg.grad, documented
% cfg.gradfile, documented
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented

% Copyright (C) 2004-2006, Robert Oostenveld
%
% $Log: dipolefitting.m,v $
% Revision 1.57  2009/07/02 15:55:15  roboos
% fixed typo
%
% Revision 1.56  2009/07/02 15:54:15  roboos
% convert 1x6 dipole position (after scanning with symmetric pair) into 2x3 to prevent warning later
%
% Revision 1.55  2009/07/02 15:37:04  roboos
% use senstype instead of strcmp
% use fixdipole helper function for consistent dipole structure representation
%
% Revision 1.54  2009/06/03 09:51:16  roboos
% changed input specification of dip.mom into 3xNdipoles
% give explicit error if unsupported dipole model
%
% Revision 1.53  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.52  2009/01/16 17:21:20  sashae
% added config tracking
%
% Revision 1.51  2008/10/02 15:32:20  sashae
% replaced call to createsubcfg with checkconfig
%
% Revision 1.50  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.49  2008/07/15 19:56:44  roboos
% moved cfg details for dipole grid to subcfg (cfg.grid)subcfg (cfg.grid.xxx)
%
% Revision 1.48  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.47  2008/03/18 12:38:40  roboos
% be more explicit about absense of symmetry constraint (since other options can also be passed in cfg.dipfit)
%
% Revision 1.46  2007/12/04 17:33:50  roboos
% added a space
%
% Revision 1.45  2007/12/04 17:33:02  roboos
% changed some print commands
%
% Revision 1.44  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.43  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.42  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.41  2007/03/21 11:25:45  roboos
% removed model=fixed, since it was largely unsupported anyway
% more consistent handling of rotating/moving dipole model ("case model" everywhere)
% more consistent handling of input dipole model specification
% improved handling of dip(t) for moving dipole model when t>1
% fixed bug as reported by Arno in the structure assignment for moving dipoles ("Subscripted assignment between dissimilar structures")
% added cfg.dipfit options to the documentation
%
% Revision 1.40  2007/01/09 09:51:17  roboos
% display the error message that was caught in try statement
%
% Revision 1.39  2006/11/23 10:52:46  roboos
% updated documentation
%
% Revision 1.38  2006/10/12 11:35:13  roboos
% replaced some ifs with switches
% reshape dipole position from 1x6 to 2x3 etc.
%
% Revision 1.37  2006/10/12 08:45:26  roboos
% added symmetry in z-direction, moved symmetric grid creation to prepare_dipole_grid
%
% Revision 1.36  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.35  2006/07/05 10:34:52  roboos
% minor change in documentation
%
% Revision 1.34  2006/07/05 10:32:14  roboos
% minor change in documentation
%
% Revision 1.33  2006/07/05 10:23:57  roboos
% updated documentation
% removed default value for cfg.resolution (now determined in prepare_dipole_grid)
% removed default 'auto' values for xgrid/ygrid/zgrid (were similar to the default values in prepare_dipole_grid)
%
% Revision 1.32  2006/07/04 17:06:23  ingnie
% updated documentation
%
% Revision 1.31  2006/06/20 16:25:58  ingnie
% updated documentation
%
% Revision 1.30  2006/06/13 14:48:08  ingnie
% updated documentation
%
% Revision 1.29  2006/05/11 07:13:51  roboos
% removed the message about which optimization function is being used
%
% Revision 1.28  2006/05/10 16:01:55  roboos
% Fixed bug that occurred when the channel ordering in the elec was
% different than the data. Added measured and model data to the output
% structure, added elec or grad to the output structure (usefull for
% topoplotting).
%
% Revision 1.27  2006/05/10 15:41:14  roboos
% Changed the order of the input arguments for the dipole_fit() function
% to be consistent with the function itself. Changed the handling of a
% failed dipole fit. Added a transparent way for handing over optional
% key-value parameters to the low-level dipole_fit() function, currently
% supported are constr, maxiter, display.
%
% Revision 1.26  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.25  2006/04/10 16:33:46  ingnie
% updated documentation
%
% Revision 1.24  2006/04/06 16:17:09  ingnie
% updated documentation
%
% Revision 1.23  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.22  2005/12/21 09:26:29  roboos
% implemented support for a fixed dipole model, only dipole moment is computed (i.e. no gridsearch or nonlinear fit is done)
% some changes in whitespace by auto indentation
%
% Revision 1.21  2005/11/01 09:37:46  roboos
% removed local subfunction compute_dipole_fit() and replaced by a
% call to the new and general dipole_fit() function which supports
% both EEG and MEG
%
% Revision 1.20  2005/10/14 16:23:11  roboos
% added support for dipole fitting in the frequency domain (based on fourier or csd input)
% restructured the handling of comp and freq input, now both consistent with timelock input
% changed the channel and latency selection of the topographies
%
% Revision 1.19  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.18  2005/06/01 09:15:36  roboos
% replaced fprintf for each grid point with more controlled progress indicator, using cfg.feedback
%
% Revision 1.17  2005/06/01 08:00:59  roboos
% some small changes to the help documentation
%
% Revision 1.16  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.15  2005/03/14 10:44:54  roboos
% added informative fprintf about optimization toolbox (no functional change)
%
% Revision 1.14  2005/02/21 07:56:54  roboos
% deleted the subfunction compute_leadfield, it is now a separate function
%
% Revision 1.13  2004/12/08 18:00:14  roboos
% implemented consistent method of selecting a subset of channels for
% forward and inverse computations using cfg.channel and updated the
% ducumentation
%
% Revision 1.12  2004/09/06 07:47:03  roboos
% fixed bug for finding optimium after grid scanning fixed dipole model
%
% Revision 1.11  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details
%
% Revision 1.10  2004/08/30 13:18:06  roboos
% fixed 2 bogs for moving dipole model
%
% Revision 1.9  2004/06/21 20:13:25  roberto
% fixed two bugs: one with sens.pnt and the other with inside/outside dipole positions after scanning
%
% Revision 1.8  2004/06/03 15:47:21  roberto
% improved detecton of correct configuration (numdipoles)
%
% Revision 1.7  2004/05/19 15:39:35  roberto
% fixed bug in time, which can be cell-array in case of ICA on raw data
% added complete version details to output configuration
%
% Revision 1.6  2004/05/13 23:35:47  roberto
% fixed multiple bugs related to new implementation moving/regional
%
% Revision 1.5  2004/05/13 22:36:30  roberto
% implemented cfg.method = 'regional' (old behaviour) and 'moving' (additional)
%
% Revision 1.4  2004/03/06 13:04:03  roberto
% fixed some small bugs and added time/component number to the output
%
% Revision 1.3  2004/02/05 12:20:06  roberto
% moved the code that was common with dipolefitting into separate prepare_grid and prepare_sens_vol functions
% fixed multiple small bugs
%
% Revision 1.2  2004/01/27 14:39:05  roberto
% multiple small bug fixes
%
% Revision 1.1  2004/01/22 21:42:52  roberto
% initial version, not tested
% support for single dipole scanning, fitting and symmetric dipole pair
%

fieldtripdefs
cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'timelock', 'freq', 'comp'}, 'feedback', 'yes');

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
cfg = checkconfig(cfg, 'createsubcfg',  {'grid'});

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
    % this structure is passed onto the low-level dipole_fit function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 -1 1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 -x1 y1 z1]
  elseif strcmp(cfg.symmetry, 'y')
    % this structure is passed onto the low-level dipole_fit function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 1 -1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 -y1 z1]
  elseif strcmp(cfg.symmetry, 'z')
    % this structure is passed onto the low-level dipole_fit function
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
if senstype(sens, 'eeg')
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

  % construct the grid on which the scanning will be done
  [grid, cfg] = prepare_dipole_grid(cfg, vol, sens);
  progress('init', cfg.feedback, 'scanning grid');
  for i=1:length(grid.inside)
    progress(i/length(grid.inside), 'scanning grid location %d/%d\n', i, length(grid.inside));
    indx = grid.inside(i);
    if isfield(grid, 'leadfield')
      % reuse the previously computed leadfield
      lf = grid.leadfield{indx};
    else
      lf = compute_leadfield(grid.pos(indx,:), sens, vol);
    end
    % the model is V=lf*mom+noise, therefore mom=pinv(lf)*V estimates the
    % dipole moment this makes the model potential U=lf*pinv(lf)*V and the
    % model error is norm(V-U) = norm(V-lf*pinv(lf)*V) = norm((eye-lf*pinv(lf))*V)
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
  progress('close');

  switch cfg.model
    case 'regional'
      % find the grid point(s) with the minimum error
      [err, indx] = min(grid.error(grid.inside));
      dip.pos = grid.pos(grid.inside(indx),:);          % note that for a symmetric dipole pair this results in a vector
      dip.pos = reshape(dip.pos, cfg.numdipoles, 3);    % convert to a Nx3 array
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
        dip(t).pos = reshape(dip(t).pos, cfg.numdipoles, 3);  % convert to a Nx3 array
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
  optarg = cfg2keyval(getfield(cfg, 'dipfit'));
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
      % instead of using dip(t) =  dipole_fit(dip(t),...), I am using temporary variables dipin and dipout
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
      lf = compute_leadfield(dip.pos, sens, vol);
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
        lf = compute_leadfield(dip(t).pos, sens, vol);
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
      % although this is technically possible sofar, it does not make any sense
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

if isfield(data, 'grad')
  % copy the gradiometer array along
  source.grad = data.grad;
end
if isfield(data, 'elec')
  % copy the electrode array along
  source.elec = data.elec;
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: dipolefitting.m,v 1.57 2009/07/02 15:55:15 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
source.cfg = cfg;

