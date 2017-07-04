function [grid, cfg] = ft_prepare_leadfield(cfg, data)

% FT_PREPARE_LEADFIELD computes the forward model for many dipole locations
% on a regular 2D or 3D grid and stores it for efficient inverse modelling
%
% Use as
%   [grid] = ft_prepare_leadfield(cfg, data);
%
% It is neccessary to input the data on which you want to perform the
% inverse computations, since that data generally contain the gradiometer
% information and information about the channels that should be included in
% the forward model computation. The data structure can be either obtained
% from FT_PREPROCESSING, FT_FREQANALYSIS or FT_TIMELOCKANALYSIS. If the data is empty,
% all channels will be included in the forward model.
%
% The configuration should contain
%   cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'),
%                            see FT_CHANNELSELECTION for details
%
% The positions of the sources can be specified as a regular 3-D
% grid that is aligned with the axes of the head coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
% Alternatively the position of a few sources at locations of interest can
% be specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = N*3 matrix with position of each source
%   cfg.grid.inside     = N*1 vector with boolean value whether grid point is inside brain (optional)
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%
% The volume conduction model of the head should be specified as
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.grad          = structure with gradiometer definition, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%   cfg.gradfile      = name of file containing the gradiometer definition, see FT_READ_SENS
%
% Optionally, you can modify the leadfields by reducing the rank (i.e.
% remove the weakest orientation), or by normalizing each column.
%   cfg.reducerank      = 'no', or number (default = 3 for EEG, 2 for MEG)
%   cfg.normalize       = 'yes' or 'no' (default = 'no')
%   cfg.normalizeparam  = depth normalization parameter (default = 0.5)
%   cfg.backproject     = 'yes' or 'no' (default = 'yes') determines when reducerank is applied
%                         whether the lower rank leadfield is projected back onto the original
%                         linear subspace, or not.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEANALYSIS, FT_DIPOLEFITTING, FT_PREPARE_HEADMODEL,
% FT_PREPARE_SOURCEMODEL

% Undocumented local options:
% cfg.feedback
% cfg.sel50p      = 'no' (default) or 'yes'
% cfg.lbex        = 'no' (default) or a number that corresponds with the radius
% cfg.mollify     = 'no' (default) or a number that corresponds with the FWHM

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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the data can be passed as input arguments or can be read from disk
hasdata = exist('data', 'var');

if ~hasdata
  % the data variable will be passed to the prepare_headmodel function below
  % where it would be used for channel selection
  data = [];
else
  % check if the input data is valid for this function
  data = ft_checkdata(data);
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});

% set the defaults
cfg.normalize      = ft_getopt(cfg, 'normalize',      'no');
cfg.normalizeparam = ft_getopt(cfg, 'normalizeparam', 0.5);
cfg.lbex           = ft_getopt(cfg, 'lbex',           'no');
cfg.sel50p         = ft_getopt(cfg, 'sel50p',         'no');
cfg.feedback       = ft_getopt(cfg, 'feedback',       'text');
cfg.mollify        = ft_getopt(cfg, 'mollify',        'no');
cfg.patchsvd       = ft_getopt(cfg, 'patchsvd',       'no');
cfg.backproject    = ft_getopt(cfg, 'backproject',    'yes'); % determines whether after rank reduction the subspace projected leadfield is backprojected onto the original space
% cfg.reducerank   = ft_getopt(cfg, 'reducerank', 'no');      % the default for this depends on EEG/MEG and is set below

% put the low-level options pertaining to the dipole grid in their own field
cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'});  % this is moved to cfg.grid.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.grid.unit by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'grid'});

% this code expects the inside to be represented as a logical array
cfg.grid = ft_checkconfig(cfg.grid, 'renamed',  {'pnt' 'pos'});
cfg = ft_checkconfig(cfg, 'index2logical', 'yes');

if strcmp(cfg.sel50p, 'yes') && strcmp(cfg.lbex, 'yes')
  ft_error('subspace projection with either lbex or sel50p is mutually exclusive');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% collect and preprocess the electrodes/gradiometer and head model
[headmodel, sens, cfg] = prepare_headmodel(cfg, data);

% set the default for reducing the rank of the leadfields
if ft_senstype(sens, 'eeg')
  cfg.reducerank = ft_getopt(cfg, 'reducerank', 3);
else
  cfg.reducerank = ft_getopt(cfg, 'reducerank', 2);
end

% construct the dipole grid according to the configuration
tmpcfg           = keepfields(cfg, {'grid', 'mri', 'headshape', 'symmetry', 'smooth', 'threshold', 'spheremesh', 'inwardshift', 'showcallinfo'});
tmpcfg.headmodel = headmodel;
tmpcfg.grad      = sens; % either electrodes or gradiometers
grid = ft_prepare_sourcemodel(tmpcfg);

% check whether units are equal (NOTE: this was previously not required,
% this check can be removed if the underlying bug is resolved. See
% http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2387
if ~isfield(headmodel, 'unit') || ~isfield(grid, 'unit') || ~isfield(sens, 'unit')
  ft_warning('cannot determine the units of all geometric objects required for leadfield computation (headmodel, sourcemodel, sensor configuration). THIS CAN LEAD TO WRONG RESULTS! (refer to http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2387)');
else
  if ~strcmp(headmodel.unit, grid.unit) || ~strcmp(grid.unit, sens.unit)
    ft_error('geometric objects (headmodel, sourcemodel, sensor configuration) are not expressed in the same units (this used to be allowed, and will be again in the future, but for now there is a bug which prevents a correct leadfield from being computed; see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2387)');
  end
end

if ft_voltype(headmodel, 'openmeeg')
  % repeated system calls to the openmeeg executable makes it rather slow
  % calling it once is much more efficient
  fprintf('calculating leadfield for all positions at once, this may take a while...\n');

  % find the indices of all grid points that are inside the brain
  insideindx = find(grid.inside);
  ndip       = length(insideindx);
  ok         = false(1,ndip);
  batchsize  = ndip;

  while ~all(ok)
    % find the first one that is not yet done
    begdip = find(~ok, 1);
    % define a batch of dipoles to jointly deal with
    enddip = min((begdip+batchsize-1), ndip); % don't go beyond the end
    batch  = begdip:enddip;
    try
      lf = ft_compute_leadfield(grid.pos(insideindx(batch),:), sens, headmodel, 'reducerank', cfg.reducerank, 'normalize', cfg.normalize, 'normalizeparam', cfg.normalizeparam);
      ok(batch) = true;
    catch
      ok(batch) = false;
      % the "catch me" syntax is broken on MATLAB74, this fixes it
      me = lasterror;
      if ~isempty(findstr(me.message, 'Output argument "dsm" (and maybe others) not assigned during call to'))
        % it does not fit in memory, split the problem in two halves and try once more
        batchsize = floor(batchsize/500);
        continue
      else
        rethrow(me);
      end % handling this particular error
    end

    % reassign the large leadfield matrix over the single grid locations
    for i=1:length(batch)
      sel = (3*i-2):(3*i);           % 1:3, 4:6, ...
      dipindx = insideindx(batch(i));
      grid.leadfield{dipindx} = lf(:,sel);
    end

    clear lf

  end % while

else
  % find the indices of all grid points that are inside the brain
  insideindx = find(grid.inside);

  ft_progress('init', cfg.feedback, 'computing leadfield');
  for i=1:length(insideindx)
    % compute the leadfield on all grid positions inside the brain
    ft_progress(i/length(insideindx), 'computing leadfield %d/%d\n', i, length(insideindx));
    thisindx = insideindx(i);
    grid.leadfield{thisindx} = ft_compute_leadfield(grid.pos(thisindx,:), sens, headmodel, 'reducerank', cfg.reducerank, 'normalize', cfg.normalize, 'normalizeparam', cfg.normalizeparam, 'backproject', cfg.backproject);

    if isfield(cfg, 'grid') && isfield(cfg.grid, 'mom')
      % multiply with the normalized dipole moment to get the leadfield in the desired orientation
      grid.leadfield{thisindx} = grid.leadfield{thisindx} * grid.mom(:,thisindx);
    end
  end % for all grid locations inside the brain
  ft_progress('close');
end

% represent the leadfield for positions outside the brain as empty array
grid.leadfield(~grid.inside) = {[]};

% add the label of the channels
grid.label           = sens.label;
grid.leadfielddimord = '{pos}_chan_ori';

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

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance grid
ft_postamble history    grid
