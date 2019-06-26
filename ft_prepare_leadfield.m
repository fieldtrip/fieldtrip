function [sourcemodel, cfg] = ft_prepare_leadfield(cfg, data)

% FT_PREPARE_LEADFIELD computes the forward model for many dipole locations
% on a regular 2D or 3D sourcemodel and stores it for efficient inverse modelling
%
% Use as
%   [sourcemodel] = ft_prepare_leadfield(cfg, data)
%
% It is necessary to input the data on which you want to perform the
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
% sourcemodel that is aligned with the axes of the head coordinate system
%   cfg.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.resolution = number (e.g. 1 cm) for automatic sourcemodel generation
% Alternatively the position of a few sources at locations of interest can
% be specified, for example obtained from an anatomical or functional MRI
%   cfg.sourcemodel.pos        = N*3 matrix with position of each source
%   cfg.sourcemodel.inside     = N*1 vector with boolean value whether sourcemodel point is inside brain (optional)
%   cfg.sourcemodel.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D sourcemodel (optional)
%
% The volume conduction model of the head should be specified as
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition or filename, see FT_READ_SENS
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
% Depending on the type of headmodel, some additional options may be
% specified.
%
% For OPENMEEG based headmodels:
%   cfg.openmeeg.batchsize    = scalar (default 100e3), number of dipoles
%                               for which the leadfield is computed in a
%                               single call to the low-level code. Trades off
%                               memory efficiency for speed.
%   cfg.openmeeg.dsm          = 'no'/'yes', reuse existing DSM if provided
%   cfg.openmeeg.keepdsm      = 'no'/'yes', option to retain DSM (no by default)
%   cfg.openmeeg.nonadaptive  = 'no'/'yes'
%
% For SINGLESHELL based headmodels:
%   cfg.singleshell.batchsize = scalar or 'all' (default 1), number of dipoles
%                               for which the leadfield is computed in a
%                               single call to the low-level code. Trades off
%                               memory efficiency for speed.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEANALYSIS, FT_DIPOLEFITTING, FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

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
cfg = ft_checkconfig(cfg, 'renamed', {'grid',    'sourcemodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'om',      'openmeeg'});
cfg = ft_checkconfig(cfg, 'renamed', {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed', {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed', {'optofile', 'opto'});

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

cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'});  % this is moved to cfg.sourcemodel.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.sourcemodel.unit by the subsequent createsubcfg

% put the low-level options pertaining to the sourcemodel in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg', {'sourcemodel'});
% move some fields from cfg.sourcemodel back to the top-level configuration
cfg = ft_checkconfig(cfg, 'createtopcfg', {'sourcemodel'});

% this code expects the inside to be represented as a logical array
cfg.sourcemodel = ft_checkconfig(cfg.sourcemodel, 'renamed',  {'pnt' 'pos'});
cfg = ft_checkconfig(cfg, 'inside2logical', 'yes');

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

% construct the sourcemodel for which the leadfield will be computed
tmpcfg           = keepfields(cfg, {'sourcemodel', 'mri', 'headshape', 'symmetry', 'smooth', 'threshold', 'spheremesh', 'inwardshift', 'xgrid' 'ygrid', 'zgrid', 'resolution', 'tight', 'warpmni', 'template', 'showcallinfo'});
tmpcfg.headmodel = headmodel;
if ft_senstype(sens, 'eeg')
  tmpcfg.elec = sens;
elseif ft_senstype(sens, 'meg')
  tmpcfg.grad = sens;
end
sourcemodel = ft_prepare_sourcemodel(tmpcfg);

% check whether units are equal (NOTE: this was previously not required,
% this check can be removed if the underlying bug is resolved. See
% http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2387
if ~isfield(headmodel, 'unit') || ~isfield(sourcemodel, 'unit') || ~isfield(sens, 'unit')
  ft_warning('cannot determine the units of all geometric objects required for leadfield computation (headmodel, sourcemodel, sensor configuration). THIS CAN LEAD TO WRONG RESULTS! (refer to http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2387)');
else
  if ~strcmp(headmodel.unit, sourcemodel.unit) || ~strcmp(sourcemodel.unit, sens.unit)
    ft_error('geometric objects (headmodel, sourcemodel, sensor configuration) are not expressed in the same units (this used to be allowed, and will be again in the future, but for now there is a bug which prevents a correct leadfield from being computed; see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2387)');
  end
end

% find the indices of all sourcemodel points that are inside the brain
insideindx = find(sourcemodel.inside);

if ft_headmodeltype(headmodel, 'openmeeg')
  
  ft_hastoolbox('openmeeg', 1);  % add to path (if not yet on path)
  
  % repeated system calls to the openmeeg executable makes it rather slow
  % calling it once is much more efficient
  fprintf('calculating leadfield for all positions at once, this may take a while...\n');

  if(~isfield(cfg,'openmeeg'))
    cfg.openmeeg = [];
  end
  batchsize   = ft_getopt(cfg.openmeeg, 'batchsize',100e3); % number of voxels per DSM batch; set to e.g. 1000 if not much RAM available
  dsm         = ft_getopt(cfg.openmeeg, 'dsm');             % reuse existing DSM if provided
  keepdsm     = ft_getopt(cfg.openmeeg, 'keepdsm', 'no');   % retain DSM
  nonadaptive = ft_getopt(cfg.openmeeg, 'nonadaptive', 'no');

  ndip       = length(insideindx);
  numchunks  = ceil(ndip/batchsize);
  if(numchunks > 1)
    if istrue(keepdsm)
      ft_warning('Keeping DSM output not supported when the computation is split into batches')
    end
    keepdsm = false;
  end

  try
    % DSM computation is computationally intensive:
    % As it can be reused with same voxel sourcemodel (i.e. if voxels are defined in
    % MRI coordinates rather than MEG coordinates), optionally save result.
    % Dense voxel grids may require several gigabytes of RAM, so optionally
    % split into smaller batches

    [h2sens,ds2sens] = ft_sensinterp_openmeeg(sourcemodel.pos(insideindx,:), headmodel, sens);

    % use pre-existing DSM if present
    if(~isempty(dsm))
      lf = ds2sens + h2sens*headmodel.mat*dsm;
    else
      lf = zeros(size(ds2sens)); % pre-allocate Msensors x Nvoxels

      for ii = 1:numchunks
        % select sourcemodel positions for this batch
        diprange = (((ii-1)*batchsize + 1):(min((ii)*batchsize,ndip)));
        % remap with 3 orientations per position
        diprangeori = [((ii-1)*3*batchsize + 1):(min((ii)*3*batchsize,3*ndip))];
        dsm = ft_sysmat_openmeeg(sourcemodel.pos(insideindx(diprange),:), headmodel, sens, nonadaptive);
        lf(:,diprangeori) = ds2sens(:,diprangeori) + h2sens*headmodel.mat*dsm;

        if istrue(keepdsm)
          % retain DSM in cfg if desired
          cfg.openmeeg.dsm = dsm;
        end

        dipindx = insideindx(diprange);
      end
    end
  catch
    me = lasterror;
    rethrow(me);
  end

  % apply montage, if applicable
  if isfield(sens, 'tra')
    lf = sens.tra * lf;
  end

  % lead field computation already done, but pass to ft_compute_leadfield so that
  % any post-computation options can be applied (e.g., normalization, etc.)
  lf = ft_compute_leadfield(sourcemodel.pos(diprange,:), sens, headmodel, 'lf', lf, 'reducerank', cfg.reducerank, 'normalize', cfg.normalize, 'normalizeparam', cfg.normalizeparam, 'backproject', cfg.backproject);

  % reshape result into sourcemodel.leadfield cell-array
  for i=1:ndip
    sourcemodel.leadfield{insideindx(i)} = lf(:,3*(i-1) + [1:3]);
  end
  clear lf

elseif ft_headmodeltype(headmodel, 'singleshell')
  cfg.singleshell = ft_getopt(cfg, 'singleshell', []);
  batchsize       = ft_getopt(cfg.singleshell, 'batchsize', 1);
  if ischar(batchsize) && strcmp(batchsize, 'all')
    batchsize = length(insideindx);
  end

  dippos     = sourcemodel.pos(insideindx,:);
  ndip       = length(insideindx);
  numchunks  = ceil(ndip/batchsize);

  ft_progress('init', cfg.feedback, 'computing leadfield');
  for k = 1:numchunks
    ft_progress(k/numchunks, 'computing leadfield %d/%d\n', k, numchunks);
    diprange = (((k-1)*batchsize + 1):(min(k*batchsize,ndip)));
    tmp      = ft_compute_leadfield(dippos(diprange,:), sens, headmodel, 'reducerank', cfg.reducerank, 'normalize', cfg.normalize, 'normalizeparam', cfg.normalizeparam, 'backproject', cfg.backproject);
    for i=1:length(diprange)
      thisindx = insideindx(diprange(i));
      if istrue(cfg.backproject)
        sourcemodel.leadfield{thisindx} = tmp(:,(i-1)*3+(1:3));
      else
        sourcemodel.leadfield{thisindx} = tmp(:,(i-1)*cfg.reducerank+(1:cfg.reducerank));
      end

      if isfield(cfg, 'sourcemodel') && isfield(cfg.sourcemodel, 'mom')
        % multiply with the normalized dipole moment to get the leadfield in the desired orientation
        sourcemodel.leadfield{thisindx} = sourcemodel.leadfield{thisindx} * sourcemodel.mom(:,thisindx);
      end
    end
  end
  ft_progress('close');

else
  ft_progress('init', cfg.feedback, 'computing leadfield');
  for i=1:length(insideindx)
    % compute the leadfield on all sourcemodel positions inside the brain
    ft_progress(i/length(insideindx), 'computing leadfield %d/%d\n', i, length(insideindx));
    thisindx = insideindx(i);
    sourcemodel.leadfield{thisindx} = ft_compute_leadfield(sourcemodel.pos(thisindx,:), sens, headmodel, 'reducerank', cfg.reducerank, 'normalize', cfg.normalize, 'normalizeparam', cfg.normalizeparam, 'backproject', cfg.backproject);

    if isfield(cfg, 'sourcemodel') && isfield(cfg.sourcemodel, 'mom')
      % multiply with the normalized dipole moment to get the leadfield in the desired orientation
      sourcemodel.leadfield{thisindx} = sourcemodel.leadfield{thisindx} * sourcemodel.mom(:,thisindx);
    end
  end % for all sourcemodel locations inside the brain
  ft_progress('close');
end

% represent the leadfield for positions outside the brain as empty array
sourcemodel.leadfield(~sourcemodel.inside) = {[]};

% add the label of the channels
sourcemodel.label           = sens.label;
sourcemodel.leadfielddimord = '{pos}_chan_ori';

% mollify the leadfields
if ~strcmp(cfg.mollify, 'no')
  sourcemodel = mollify(cfg, sourcemodel);
end

% combine leadfields in patches and do an SVD on them
if ~strcmp(cfg.patchsvd, 'no')
  sourcemodel = patchsvd(cfg, sourcemodel);
end

% compute the 50 percent channel selection subspace projection
if ~strcmp(cfg.sel50p, 'no')
  sourcemodel = sel50p(cfg, sourcemodel, sens);
end

% compute the local basis function expansion (LBEX) subspace projection
if ~strcmp(cfg.lbex, 'no')
  sourcemodel = lbex(cfg, sourcemodel);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance sourcemodel
ft_postamble history    sourcemodel
ft_postamble savevar    sourcemodel
