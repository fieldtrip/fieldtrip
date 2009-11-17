function [grid, cfg] = prepare_leadfield(cfg, data)

% PREPARE_LEADFIELD computes the forward model for many dipole locations
% on a regular 2D or 3D grid and stores it for efficient inverse modelling
%
% Use as
%   [grid] = prepare_leadfield(cfg, data);
%
% It is neccessary to input the data on which you want to perform the
% inverse computations, since that data generally contain the gradiometer
% information and information about the channels that should be included in
% the forward model computation. The data structure can be either obtained
% from PREPROCESSING, FREQANALYSIS or TIMELOCKANALYSIS. If the data is empty,
% all channels will be included in the forward model.
%
% The configuration should contain
%   cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'),
%                            see CHANNELSELECTION for details
%
% The positions of the sources can be specified as a regular 3-D
% grid that is aligned with the axes of the head coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
% Alternatively the position of a few sources at locations of interest can
% be specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = Nx3 matrix with position of each source
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%   cfg.grid.inside     = vector with indices of the sources inside the brain (optional)
%   cfg.grid.outside    = vector with indices of the sources outside the brain (optional)
%
% You should specify the volume conductor model with
%   cfg.hdmfile         = string, file containing the volume conduction model
% or alternatively
%   cfg.vol             = structure with volume conduction model
%
% If the sensor information is not contained in the data itself you should
% also specify the sensor information using
%   cfg.gradfile        = string, file containing the gradiometer definition
%   cfg.elecfile        = string, file containing the electrode definition
% or alternatively
%   cfg.grad            = structure with gradiometer definition
%   cfg.elec            = structure with electrode definition
%
% Optionally, you can modify the leadfields by reducing the rank (i.e.
% remove the weakest orientation), or by normalizing each column.
%   cfg.reducerank      = 'no', or number (default = 3 for EEG, 2 for MEG)
%   cfg.normalize       = 'yes' or 'no' (default = 'no')
%   cfg.normalizeparam  = depth normalization parameter (default = 0.5)
%
% See also SOURCEANALYSIS

% Undocumented local options:
% cfg.feedback
% cfg.sel50p      = 'no' (default) or 'yes'
% cfg.lbex        = 'no' (default) or a number that corresponds with the radius
% cfg.mollify     = 'no' (default) or a number that corresponds with the FWHM

% This function depends on PREPARE_DIPOLE_GRID which has the following options:
% cfg.grid.xgrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.xgrid = 'auto'), documented
% cfg.grid.ygrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.ygrid = 'auto'), documented
% cfg.grid.zgrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.zgrid = 'auto'), documented
% cfg.grid.resolution, documented
% cfg.grid.pos, documented
% cfg.grid.dim, documented
% cfg.grid.inside, documented
% cfg.grid.outside, documented
% cfg.mri
% cfg.mriunits
% cfg.smooth
% cfg.sourceunits
% cfg.threshold
%
% This function depends on PREPARE_VOL_SENS which has the following options:
% cfg.channel, documented
% cfg.elec, documented
% cfg.elecfile, documented
% cfg.grad, documented
% cfg.gradfile, documented
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented

% Copyright (C) 2004-2006, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs
cfg = checkconfig(cfg, 'trackconfig', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2
  data = [];
end

% set the defaults
if ~isfield(cfg, 'normalize'),        cfg.normalize  = 'no';          end
if ~isfield(cfg, 'normalizeparam'),   cfg.normalizeparam = 0.5;       end
if ~isfield(cfg, 'lbex'),             cfg.lbex       = 'no';          end
if ~isfield(cfg, 'sel50p'),           cfg.sel50p     = 'no';          end
if ~isfield(cfg, 'feedback'),         cfg.feedback   = 'text';        end
if ~isfield(cfg, 'mollify'),          cfg.mollify    = 'no';          end
if ~isfield(cfg, 'patchsvd'),         cfg.patchsvd   = 'no';          end
% if ~isfield(cfg, 'reducerank'),     cfg.reducerank = 'no';          end  % the default for this depends on EEG/MEG and is set below

% put the low-level options pertaining to the dipole grid in their own field
cfg = checkconfig(cfg, 'createsubcfg',  {'grid'});

if strcmp(cfg.sel50p, 'yes') && strcmp(cfg.lbex, 'yes')
  error('subspace projection with either lbex or sel50p is mutually exclusive');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect and preprocess the electrodes/gradiometer and head model
[vol, sens, cfg] = prepare_headmodel(cfg, data);

% set the default for reducing the rank of the leadfields
if ~isfield(cfg, 'reducerank')
  if senstype(sens, 'eeg')
    cfg.reducerank = 3;
  else
    cfg.reducerank = 2;
  end
end

% construct the grid on which the scanning will be done
[grid, cfg] = prepare_dipole_grid(cfg, vol, sens);

progress('init', cfg.feedback, 'computing leadfield');
if strcmp(vol.type, 'openmeeg') 
  fprintf('Calculating infinite medium solution for all positions\n Please Wait...\n');
  dsm = openmeeg_dsm(grid.pos,vol);
  for i=1:length(grid.inside)
    % compute the leadfield on all grid positions inside the brain
    progress(i/length(grid.inside), 'computing leadfield %d/%d\n',i, length(grid.inside));    
    dipindx = grid.inside(i);
    tmp = vol.mat*dsm(:,[(3*(i-1)):(3*(i-1)+2)]+1);
    grid.leadfield{dipindx} = tmp - repmat(mean(tmp,1),size(tmp,1),1);
    if isfield(cfg, 'grid') && isfield(cfg.grid, 'mom')
      % multiply with the normalized dipole moment to get the leadfield in the desired orientation
      grid.leadfield{dipindx} = grid.leadfield{dipindx} * grid.mom(:,dipindx);
    end   
  end
  progress('close');

else
  for i=1:length(grid.inside)
    % compute the leadfield on all grid positions inside the brain
    progress(i/length(grid.inside), 'computing leadfield %d/%d\n', i, length(grid.inside));
    dipindx = grid.inside(i);
    grid.leadfield{dipindx} = compute_leadfield(grid.pos(dipindx,:), sens, vol, 'reducerank', cfg.reducerank, 'normalize', cfg.normalize, 'normalizeparam', cfg.normalizeparam);
    
    if isfield(cfg, 'grid') && isfield(cfg.grid, 'mom')
      % multiply with the normalized dipole moment to get the leadfield in the desired orientation
      grid.leadfield{dipindx} = grid.leadfield{dipindx} * grid.mom(:,dipindx);
    end
  end % for all grid locations inside the brain
  progress('close');
end


% fill the positions outside the brain with NaNs
grid.leadfield(grid.outside) = {nan};

% mollify the leadfields
if ~strcmp(cfg.mollify, 'no')
  grid = mollify(cfg, grid);
end

% combine leadfields in patches and do an SVD on them
if ~strcmp(cfg.patchsvd, 'no')
  grid = patchsvd(cfg, grid);
end

% compute the 50 percent channel selection subspace projection
if ~strcmp(cfg.sel50p, 'no')
  grid = sel50p(cfg, grid, sens);
end

% compute the local basis function expansion (LBEX) subspace projection
if ~strcmp(cfg.lbex, 'no')
  grid = lbex(cfg, grid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: prepare_leadfield.m,v 1.30 2009/04/09 16:28:40 crimic Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
grid.cfg = cfg;

