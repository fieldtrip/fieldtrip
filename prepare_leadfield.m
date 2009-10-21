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
% $Log: prepare_leadfield.m,v $
% Revision 1.30  2009/04/09 16:28:40  crimic
% updated documentation according to guidelines
%
% Revision 1.29  2009/03/23 13:41:56  roboos
% use senstype instead of sensortype
%
% Revision 1.28  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.27  2009/01/16 17:21:20  sashae
% added config tracking
%
% Revision 1.26  2008/10/02 15:32:21  sashae
% replaced call to createsubcfg with checkconfig
%
% Revision 1.25  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.24  2008/07/15 19:56:44  roboos
% moved cfg details for dipole grid to subcfg (cfg.grid)subcfg (cfg.grid.xxx)
%
% Revision 1.23  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.22  2008/03/18 21:14:14  roboos
% use grid.mom instead of grid.ori, mom should also be 3xNdips
%
% Revision 1.21  2007/11/05 09:47:44  roboos
% added cfg.normalizeparam for depth normalization
%
% Revision 1.20  2007/05/15 07:01:55  roboos
% updated help
%
% Revision 1.19  2006/10/04 08:10:51  roboos
% changed default for cfg.reducerank, new default is 2 for MEG and 3 for EEG
%
% Revision 1.18  2006/09/19 11:03:28  jansch
% added support for mollification and patchsvd
%
% Revision 1.17  2006/07/24 07:59:16  roboos
% updated documentation
%
% Revision 1.16  2006/07/05 10:32:15  roboos
% minor change in documentation
%
% Revision 1.15  2006/07/05 10:26:27  roboos
% updated documentation, the cfg.xgrid/ygrid/zgrid options have moved to cfg.grid substructure
%
% Revision 1.14  2006/05/23 10:18:27  roboos
% return the cfg as second argument
%
% Revision 1.13  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.12  2006/04/10 16:33:46  ingnie
% updated documentation
%
% Revision 1.11  2006/04/06 16:17:10  ingnie
% updated documentation
%
% Revision 1.10  2005/11/08 11:21:58  roboos
% moved the handling of the normalize and reducerank options towards compute_leadfield
%
% Revision 1.9  2005/10/14 16:23:48  roboos
% fixed typo in fprintf (leedfield->leadfield)
%
% Revision 1.8  2005/09/29 12:52:25  jansch
% changed column-wise normalisation into matrix-wise normalisation
%
% Revision 1.7  2005/09/29 10:37:08  jansch
% fixed bug in improper normalisation
%
% Revision 1.6  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.5  2005/06/08 13:40:33  roboos
% replaced the specific call to either meg_leadfield or eeg_leadfield to the generic compute_leadfield
%
% Revision 1.4  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.3  2005/04/25 16:05:11  roboos
% includedupdated options for cortex segmentation based grid in the documentation (not in main documentation)
%
% Revision 1.2  2005/01/10 16:28:10  roboos
% some changes in order of code, no functional change
%
% Revision 1.1  2004/12/24 14:31:09  roboos
% Renamed precompute_leadfield into prepare_leadfield. For backward
% compatibility in peoples scripts, precompute_leadfield still will remain
% available for some time (it prints a warning and calls prepare_leadfield).
% The newly committed prepare_leadfield is a copy of precompute_leadfield
% revision 1.6.
%
% Revision 1.6  2004/12/08 18:00:13  roboos
% implemented consistent method of selecting a subset of channels for
% forward and inverse computations using cfg.channel and updated the
% ducumentation
%
% Revision 1.5  2004/10/27 16:13:53  roboos
% changed progress indication into private progress function
% added sel50p subspace projection
%
% Revision 1.4  2004/09/22 08:49:12  roboos
% added computation of LBEX
%
% Revision 1.3  2004/08/19 06:31:31  roboos
% removed consistency check for cfg.xgrid etc., this is now done in prepare_dipole_grid
%
% Revision 1.2  2004/08/17 08:53:55  roboos
% added cfg to the output, updated help
%
% Revision 1.1  2004/08/16 09:10:47  roboos
% initial version, this replaces the leadfield precomputation in sourceanalysis with some extra options
%

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

