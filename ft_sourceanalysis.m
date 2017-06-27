function [source] = ft_sourceanalysis(cfg, data, baseline)

% FT_SOURCEANALYSIS performs beamformer dipole analysis on EEG or MEG data
% after preprocessing and a timelocked or frequency analysis
%
% Use as either
%   [source] = ft_sourceanalysis(cfg, freq)
%   [source] = ft_sourceanalysis(cfg, timelock)
%
% where the data in freq or timelock should be organised in a structure
% as obtained from the FT_FREQANALYSIS or FT_TIMELOCKANALYSIS function. The
% configuration "cfg" is a structure containing information about
% source positions and other options.
%
% The different source reconstruction algorithms that are implemented
% are
%   cfg.method     = 'lcmv'    linear constrained minimum variance beamformer
%                    'sam'     synthetic aperture magnetometry
%                    'dics'    dynamic imaging of coherent sources
%                    'pcc'     partial cannonical correlation/coherence
%                    'mne'     minimum norm estimation
%                    'rv'      scan residual variance with single dipole
%                    'music'   multiple signal classification
%                    'sloreta' standardized low-resolution electromagnetic tomography
%                    'eloreta' exact low-resolution electromagnetic tomography
% The DICS and PCC methods are for frequency or time-frequency domain data, all other
% methods are for time domain data. ELORETA can be used both for time, frequency and
% time-frequency domain data.
%
% The source model to use in the reconstruction should be specified as
%   cfg.grid            = structure, see FT_PREPARE_SOURCEMODEL or FT_PREPARE_LEADFIELD
% The positions of the dipoles can be specified as a regular 3-D
% grid that is aligned with the axes of the head coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
%   cfg.grid.inside     = N*1 vector with boolean value whether grid point is inside brain (optional)
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
% If the source model destribes a triangulated cortical sheet, it is described as
%   cfg.grid.pos        = N*3 matrix with the vertex positions of the cortical sheet
%   cfg.grid.tri        = M*3 matrix that describes the triangles connecting the vertices
% Alternatively the position of a few dipoles at locations of interest can be
% specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = N*3 matrix with position of each source
%
% Besides the source positions, you may also include previously computed
% spatial filters and/or leadfields like this
%   cfg.grid.filter
%   cfg.grid.leadfield
%
% The following strategies are supported to obtain statistics for the source parameters using
% multiple trials in the data, either directly or through a resampling-based approach
%   cfg.rawtrial      = 'no' or 'yes'   construct filter from single trials, apply to single trials. Note that you also may want to set cfg.keeptrials='yes' to keep all trial information, especially if using in combination with grid.filter
%   cfg.jackknife     = 'no' or 'yes'   jackknife resampling of trials
%   cfg.pseudovalue   = 'no' or 'yes'   pseudovalue resampling of trials
%   cfg.bootstrap     = 'no' or 'yes'   bootstrap resampling of trials
%   cfg.numbootstrap  = number of bootstrap replications (e.g. number of original trials)
% If none of these options is specified, the average over the trials will
% be computed prior to computing the source reconstruction.
%
% To obtain statistics over the source parameters between two conditions, you
% can also use a resampling procedure that reshuffles the trials over both
% conditions. In that case, you should call the function with two datasets
% containing single trial data like
%   [source] = ft_sourceanalysis(cfg, freqA, freqB)
%   [source] = ft_sourceanalysis(cfg, timelockA, timelockB)
% and you should specify
%   cfg.randomization      = 'no' or 'yes'
%   cfg.permutation        = 'no' or 'yes'
%   cfg.numrandomization   = number, e.g. 500
%   cfg.numpermutation     = number, e.g. 500 or 'all'
%
% If you have not specified a grid with pre-computed leadfields,
% the leadfield for each grid location will be computed on the fly.
% In that case you can modify the leadfields by reducing the rank
% (i.e.  remove the weakest orientation), or by normalizing each
% column.
%   cfg.reducerank  = 'no', or number (default = 3 for EEG, 2 for MEG)
%   cfg.normalize   = 'no' or 'yes' (default = 'no')
%
% Other configuration options are
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see FT_CHANNELSELECTION for details
%   cfg.frequency     = single number (in Hz)
%   cfg.latency       = single number in seconds, for time-frequency analysis
%   cfg.lambda        = number or empty for automatic default
%   cfg.refchan       = reference channel label (for coherence)
%   cfg.refdip        = reference dipole location (for coherence)
%   cfg.supchan       = suppressed channel label(s)
%   cfg.supdip        = suppressed dipole location(s)
%   cfg.keeptrials    = 'no' or 'yes'
%   cfg.keepleadfield = 'no' or 'yes'
%   cfg.projectnoise  = 'no' or 'yes'
%   cfg.keepfilter    = 'no' or 'yes'
%   cfg.keepcsd       = 'no' or 'yes'
%   cfg.keepmom       = 'no' or 'yes'
%   cfg.feedback      = 'no', 'text', 'textbar', 'gui' (default = 'text')
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
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_SOURCEDESCRIPTIVES, FT_SOURCESTATISTICS, FT_PREPARE_LEADFIELD,
% FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

% Undocumented local options:
% cfg.numcomponents
% cfg.refchannel
% cfg.trialweight   = 'equal' or 'proportional'
% cfg.powmethod     = 'lambda1' or 'trace'

% Copyright (c) 2003-2008, Robert Oostenveld, F.C. Donders Centre
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
ft_preamble loadvar data baseline
ft_preamble provenance data baseline
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the baseline data can be passed as input argument or can be read from disk
hasbaseline = exist('baseline', 'var');

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'timelock', 'freq', 'comp'}, 'feedback', 'yes');

if hasbaseline
  baseline = ft_checkdata(baseline, 'datatype', {'timelock', 'freq', 'comp'}, 'feedback', 'yes');
end

% check that the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',     {'toilim', 'latency'});
cfg = ft_checkconfig(cfg, 'renamed',     {'foilim', 'frequency'});
cfg = ft_checkconfig(cfg, 'renamed',     {'jacknife',   'jackknife'});
cfg = ft_checkconfig(cfg, 'renamed',     {'refchannel', 'refchan'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'power',           'dics'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'coh_refchan',     'dics'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'coh_refdip',      'dics'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'dics_cohrefchan', 'dics'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'dics_cohrefdip',  'dics'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'parallel', 'trials'});
cfg = ft_checkconfig(cfg, 'forbidden', {'foi', 'toi'});
cfg = ft_checkconfig(cfg, 'renamed',     {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed',     {'vol',     'headmodel'});

% determine the type of input data
isfreq     = ft_datatype(data, 'freq');
iscomp     = ft_datatype(data, 'comp');
istimelock = ft_datatype(data, 'timelock');
if all(~[isfreq iscomp istimelock])
  error('input data is not recognized');
end

% set the defaults
cfg.keeptrials       = ft_getopt(cfg, 'keeptrials', 'no');
cfg.keepleadfield    = ft_getopt(cfg, 'keepleadfield', 'no');
cfg.trialweight      = ft_getopt(cfg, 'trialweight', 'equal');
cfg.jackknife        = ft_getopt(cfg, 'jackknife',   'no');
cfg.pseudovalue      = ft_getopt(cfg, 'pseudovalue', 'no');
cfg.bootstrap        = ft_getopt(cfg, 'bootstrap',   'no');
cfg.singletrial      = ft_getopt(cfg, 'singletrial', 'no');
cfg.rawtrial         = ft_getopt(cfg, 'rawtrial',    'no');
cfg.randomization    = ft_getopt(cfg, 'randomization', 'no');
cfg.numrandomization = ft_getopt(cfg, 'numrandomization', 100);
cfg.permutation      = ft_getopt(cfg, 'permutation',      'no');
cfg.numpermutation   = ft_getopt(cfg, 'numpermutation',   100);
cfg.wakewulf         = ft_getopt(cfg, 'wakewulf', 'yes');
cfg.killwulf         = ft_getopt(cfg, 'killwulf', 'yes');
cfg.channel          = ft_getopt(cfg, 'channel',  'all');
cfg.supdip           = ft_getopt(cfg, 'supdip',        []);
cfg.latency          = ft_getopt(cfg, 'latency',   'all');
cfg.frequency        = ft_getopt(cfg, 'frequency', 'all');

if istimelock
  cfg.method = ft_getopt(cfg, 'method', 'lcmv');
elseif isfreq
  cfg.method = ft_getopt(cfg, 'method', 'dics');
else
  cfg.method = ft_getopt(cfg, 'method', []);
end

% put the low-level options pertaining to the source reconstruction method in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg',  cfg.method);

% put the low-level options pertaining to the dipole grid in their own field
cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'});  % this is moved to cfg.grid.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.grid.unit  by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'createsubcfg', 'grid');

cfg.(cfg.method).keepfilter    = ft_getopt(cfg.(cfg.method), 'keepfilter',    'no');
cfg.(cfg.method).keepcsd       = ft_getopt(cfg.(cfg.method), 'keepcsd',       'no');
cfg.(cfg.method).keepmom       = ft_getopt(cfg.(cfg.method), 'keepmom',       'yes');
cfg.(cfg.method).projectnoise  = ft_getopt(cfg.(cfg.method), 'projectnoise',  'no');
cfg.(cfg.method).feedback      = ft_getopt(cfg.(cfg.method), 'feedback',      'text');
cfg.(cfg.method).lambda        = ft_getopt(cfg.(cfg.method), 'lambda',        []);
cfg.(cfg.method).powmethod     = ft_getopt(cfg.(cfg.method), 'powmethod',     []);
cfg.(cfg.method).normalize     = ft_getopt(cfg.(cfg.method), 'normalize',     'no');
cfg.(cfg.method).reducerank     = ft_getopt(cfg.(cfg.method), 'reducerank',    []); % the default for this is handled below

if hasbaseline && (strcmp(cfg.randomization, 'no') && strcmp(cfg.permutation, 'no'))
  error('input of two conditions only makes sense if you want to randomize or permute');
elseif ~hasbaseline && (strcmp(cfg.randomization, 'yes') || strcmp(cfg.permutation, 'yes'))
  error('randomization or permutation requires that you give two conditions as input');
end

if isfield(cfg, 'latency') && ischar(cfg.latency) && strcmp(cfg.latency, 'all') && istimelock
  %error('specification of cfg.latency is only required for time-frequency data');
end

if sum([strcmp(cfg.jackknife, 'yes'), strcmp(cfg.bootstrap, 'yes'), strcmp(cfg.pseudovalue, 'yes'), strcmp(cfg.singletrial, 'yes'), strcmp(cfg.rawtrial, 'yes'), strcmp(cfg.randomization, 'yes'), strcmp(cfg.permutation, 'yes')])>1
  error('jackknife, bootstrap, pseudovalue, singletrial, rawtrial, randomization and permutation are mutually exclusive');
end

if strcmp(cfg.rawtrial,'yes') && isfield(cfg,'grid') && ~isfield(cfg.grid,'filter')
  warning('Using each trial to compute its own filter is not currently recommended. Use this option only with precomputed filters in grid.filter');
end

% start with an empty output structure
source = [];

if istimelock
  % add the time axis to the output
  tmpcfg = keepfields(cfg, {'channel', 'latency', 'showcallinfo'});
  tmpcfg.avgovertime = 'no';
  data = ft_selectdata(tmpcfg, data);
  % restore the provenance information
  [cfg, data] = rollback_provenance(cfg, data);
  
  % copy the descriptive fields to the output
  source = copyfields(data, source, {'time'});
  
elseif isfreq
  tmpcfg = keepfields(cfg, {'channel', 'latency', 'frequency', 'refchan', 'nanmean', 'showcallinfo'});
  
  % ensure that the refchan is kept, if present
  if isfield(tmpcfg, 'refchan') && ~isempty(tmpcfg.refchan) && isempty(match_str(tmpcfg.channel, tmpcfg.refchan))
    hasrefchan = 1;
  else
    hasrefchan = 0;
  end
  
  if hasrefchan,
    if ischar(tmpcfg.refchan), tmpcfg.refchan = {tmpcfg.refchan}; end
    tmpchannel     = ft_channelselection(tmpcfg.channel, data.label); % the channels needed for the spatial filter
    tmpcfg.channel = cat(1, tmpchannel, tmpcfg.refchan);
    tmpcfg         = rmfield(tmpcfg, 'refchan');
  end
  
  tmpcfg.avgoverfreq = 'yes';
  if isfield(data, 'time')
    tmpcfg.avgovertime = 'yes';
  end
  data = ft_selectdata(tmpcfg, data);
  % restore the provenance information
  [cfg, data] = rollback_provenance(cfg, data);
  
  if hasrefchan, cfg.channel = match_str(data.label, tmpchannel); end
  
  % copy the descriptive fields to the output
  source = copyfields(data, source, {'time', 'freq', 'cumtapcnt'});
  
  % HACK the remainder of the code expects a single number
  cfg.frequency = mean(cfg.frequency);
  if isfield(data, 'time')
    cfg.latency   = mean(cfg.latency);
  end
  
elseif iscomp
  % FIXME, select the components here
  % FIXME, add the component numbers to the output
  error('the use of component data in ft_sourceanalysis is disabled for the time being: if you encounter this error message and you need this functionality please contact the FieldTrip development team');
end

convertcomp = false;
if iscomp && (strcmp(cfg.method, 'rv') || strcmp(cfg.method, 'music'))
  % these timelock methods are also supported for frequency or component data
  if iscomp
    convertcomp = true;
    % the conversion will be done below, after the latency and channel selection
  end
elseif isfreq && isfield(data, 'labelcmb')
  % ensure that the cross-spectral densities are chan_chan_therest,
  % otherwise the latency and frequency selection could fail, so we don't
  % need to worry about linearly indexed cross-spectral densities below
  % this point, this step may take some time, if multiple trials are
  % present in the data
  fprintf('converting the linearly indexed channelcombinations into a square CSD-matrix\n');
  data = ft_checkdata(data, 'cmbrepresentation', 'full');
end

if isfreq
  % as per the call to ft_checkdata above, the dimord of the freq-data is
  % either (rpt_)chan_chan_otherstuff, or rpttap_chan_otherstuff. The
  % former is with cross-spectra, the latter is with fourierspctrm
  
  % previously there was some explicit dimord checking here, but I think
  % with the more consistent data handling it is not necessary anymore.
end

% collect and preprocess the electrodes/gradiometer and head model
[headmodel, sens, cfg] = prepare_headmodel(cfg, data);

% set the default for reducing the rank of the leadfields
if isempty(cfg.(cfg.method).reducerank)
  if ft_senstype(sens, 'eeg')
    cfg.(cfg.method).reducerank = 'no';    % for EEG
  elseif ft_senstype(sens, 'meg') && ft_voltype(headmodel, 'infinite')
    cfg.(cfg.method).reducerank = 'no';    % for MEG with a magnetic dipole, e.g. a HPI coil
  elseif ft_senstype(sens, 'meg')
    cfg.(cfg.method).reducerank = 'yes';   % for MEG with a current dipole in a volume conductor
  end
end

% It might be that the number of channels in the data, the number of
% channels in the electrode/gradiometer definition and the number of
% channels in the localspheres volume conduction model are different.
% Hence a subset of the data channels will be used.
Nchans = length(cfg.channel);

if strcmp(cfg.keepleadfield, 'yes') && (~isfield(cfg, 'grid') || ~isfield(cfg.grid, 'leadfield'))
  % precompute the leadfields upon the users request
  fprintf('precomputing leadfields\n');
  grid = ft_prepare_leadfield(cfg, data);
elseif (strcmp(cfg.permutation,   'yes') || ...
    strcmp(cfg.randomization, 'yes') || ...
    strcmp(cfg.bootstrap,     'yes') || ...
    strcmp(cfg.jackknife,      'yes') || ...
    strcmp(cfg.pseudovalue,   'yes') || ...
    strcmp(cfg.singletrial,   'yes') || ...
    strcmp(cfg.rawtrial,      'yes')) && (~isfield(cfg, 'grid') || ~isfield(cfg.grid, 'leadfield'))
  % also precompute the leadfields if multiple trials have to be processed
  fprintf('precomputing leadfields for efficient handling of multiple trials\n');
  grid = ft_prepare_leadfield(cfg, data);
else
  % only prepare the dipole grid positions, the leadfield will be computed on the fly if not present
  
  % copy all options that are potentially used in ft_prepare_sourcemodel
  tmpcfg           = keepfields(cfg, {'grid' 'mri' 'headshape' 'symmetry' 'smooth' 'threshold' 'spheremesh' 'inwardshift', 'showcallinfo'});
  tmpcfg.headmodel = headmodel;
  tmpcfg.grad      = sens; % this can be electrodes or gradiometers
  grid = ft_prepare_sourcemodel(tmpcfg);
  
end

if isfield(cfg.grid, 'filter')
  if numel(cfg.grid.filter) == size(grid.pos, 1)
    grid.filter = cfg.grid.filter;
  else
    ft_warning('ignoring predefined filter as it does not match the number of source positions');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do frequency domain source reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfreq && any(strcmp(cfg.method, {'dics', 'pcc', 'eloreta', 'mne','harmony', 'rv', 'music'}))
  
  switch cfg.method
    case 'pcc'
      % this can handle both a csd matrix and a fourier matrix, but needs
      % special handling of the refdip etc.
      cfg.refdip = ft_getopt(cfg, 'refdip', []);
      cfg.supdip = ft_getopt(cfg, 'supdip', []);
      cfg.refchan = ft_getopt(cfg, 'refchan', []);
      cfg.supchan = ft_getopt(cfg, 'supchan', []);
      cfg.refchan = ft_channelselection(cfg.refchan, data.label);
      cfg.supchan = ft_channelselection(cfg.supchan, data.label);
      
      % HACK: use some experimental code
      if hasbaseline
        error('not supported')
      end
      
      tmpcfg         = cfg;
      tmpcfg.refchan = ''; % prepare_freq_matrices should not know explicitly about the refchan
      tmpcfg.channel = cfg.channel(:)';
      if isfield(cfg, 'refchan')
        % add the refchan implicitely
        tmpcfg.channel = [tmpcfg.channel cfg.refchan(:)'];
      end
      if isfield(cfg, 'supchan')
        % add the supchan implicitely
        tmpcfg.channel = [tmpcfg.channel cfg.supchan(:)'];
      end
      
      % select the data in the channels and the frequency of interest
      [Cf, Cr, Pr, Ntrials, tmpcfg] = prepare_freq_matrices(tmpcfg, data);
      
      if isfield(cfg, 'refchan') && ~isempty(cfg.refchan)
        [dum, refchanindx] = match_str(cfg.refchan, tmpcfg.channel);
      else
        refchanindx = [];
      end
      if isfield(cfg, 'supchan') && ~isempty(cfg.supchan)
        [dum, supchanindx] = match_str(cfg.supchan,tmpcfg.channel);
      else
        supchanindx = [];
      end
      Nchans = length(tmpcfg.channel); % update the number of channels
      
      % if the input data has a complete fourier spectrum, it can be projected through the filters
      if isfield(data, 'fourierspctrm')
        [dum, datchanindx] = match_str(tmpcfg.channel, data.label);
        fbin = nearest(data.freq, cfg.frequency);
        if strcmp(data.dimord, 'chan_freq')
          avg = data.fourierspctrm(datchanindx, fbin);
        elseif strcmp(data.dimord, 'rpt_chan_freq') || strcmp(data.dimord, 'rpttap_chan_freq'),
          avg = transpose(data.fourierspctrm(:, datchanindx, fbin));
        elseif strcmp(data.dimord, 'chan_freq_time')
          tbin = nearest(data.time, cfg.latency);
          avg = data.fourierspctrm(datchanindx, fbin, tbin);
        elseif strcmp(data.dimord, 'rpt_chan_freq_time') || strcmp(data.dimord, 'rpttap_chan_freq_time'),
          tbin = nearest(data.time, cfg.latency);
          avg  = transpose(data.fourierspctrm(:, datchanindx, fbin, tbin));
        end
      else
        avg = [];
      end
      
    case {'eloreta' 'mne' 'rv' 'music' 'harmony'}
      % these can handle both a csd matrix and a fourier matrix
      [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices(cfg, data);
      
      % if the input data has a complete fourier spectrum, it can be projected through the filters
      if isfield(data, 'fourierspctrm')
        [dum, datchanindx] = match_str(cfg.channel, data.label);
        fbin = nearest(data.freq, cfg.frequency);
        if strcmp(data.dimord, 'chan_freq')
          avg = data.fourierspctrm(datchanindx, fbin);
        elseif strcmp(data.dimord, 'rpt_chan_freq') || strcmp(data.dimord, 'rpttap_chan_freq'),
          avg = transpose(data.fourierspctrm(:, datchanindx, fbin));
        elseif strcmp(data.dimord, 'chan_freq_time')
          tbin = nearest(data.time, cfg.latency);
          avg = data.fourierspctrm(datchanindx, fbin, tbin);
        elseif strcmp(data.dimord, 'rpt_chan_freq_time') || strcmp(data.dimord, 'rpttap_chan_freq_time'),
          tbin = nearest(data.time, cfg.latency);
          avg  = transpose(data.fourierspctrm(:, datchanindx, fbin, tbin));
        end
      else % The input data is a CSD matrix, this is enough for computing source power, coherence and residual power.
        avg = Cf;
      end
      
    case 'dics'
      
      [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices(cfg, data);
      
      % assign a descriptive name to each of the dics sub-methods, the default is power only
      if isfield(cfg, 'refdip') && ~isempty(cfg.refdip);
        submethod = 'dics_refdip';
      elseif isfield(cfg, 'refchan') && ~isempty(cfg.refchan);
        submethod = 'dics_refchan';
      else
        submethod = 'dics_power';
      end
      
    otherwise
  end
  
  % This is the place to check for the consistency of the channel order in
  % the pre-computed leadfields/spatial filters, and to correct for it, if
  % necessary. This pertains to bugs 1746 and 3029.
  if isfield(grid, 'label') && (isfield(grid, 'leadfield') || isfield(grid, 'filter'))
    % match the channels in the leadfields/filters with those in the data
    [i1, i2] = match_str(cfg.channel, grid.label);
    if ~isequal(i2(:), (1:numel(grid.label))')
      if isfield(grid, 'leadfield')
        fprintf('\n\nSubselecting/reordering the channels in the precomputed leadfields\n\n');
        inside_indx = find(grid.inside);
        for k = inside_indx(:)'
          grid.leadfield{k} = grid.leadfield{k}(i2, :);
        end
      end
      if isfield(grid, 'filter')
        fprintf('\n\nSubselecting/reordering the channels in the precomputed filters\n\n');
        inside_indx = find(grid.inside);
        for k = inside_indx(:)'
          grid.filter{k} = grid.filter{k}(:, i2);
        end
      end
      grid.label = grid.label(i2);
    end
    if ~isequal(i1(:), (1:numel(cfg.channel))')
      % this is not so easy to deal with, throw an error
      error('There''s a mismatch between the number/order of channels in the data, with respect to the channels in the precomputed leadfield/filter. This is not easy to solve automatically. Please look into this.');
    end
  end
  
  
  % fill these with NaNs, so that I dont have to treat them separately
  if isempty(Cr), Cr = nan(Ntrials, Nchans, 1); end
  if isempty(Pr), Pr = nan(Ntrials, 1, 1); end
  
  if hasbaseline
    % repeat the conversion for the baseline condition
    [bCf, bCr, bPr, Nbaseline, cfg] = prepare_freq_matrices(cfg, baseline);
    % fill these with NaNs, so that I dont have to treat them separately
    if isempty(bCr), bCr = nan(Nbaseline, Nchans, 1); end
    if isempty(bPr), bPr = nan(Nbaseline, 1, 1); end
    % rename the active condition for convenience
    aCf = Cf;
    aCr = Cr;
    aPr = Pr;
    % this is required for averaging 2 conditions using prepare_resampled_data
    cfg2 = [];
    cfg2.numcondition = 2;
    % this is required for randomizing/permuting 2 conditions using prepare_resampled_data
    cfg.numcondition = 2;
  end
  
  % prepare the resampling of the trials, or average the data if multiple trials are present and no resampling is neccessary
  if (Ntrials<=1) && (strcmp(cfg.jackknife, 'yes') || strcmp(cfg.bootstrap, 'yes') || strcmp(cfg.pseudovalue, 'yes') || strcmp(cfg.singletrial, 'yes') || strcmp(cfg.rawtrial, 'yes') || strcmp(cfg.randomization, 'yes') || strcmp(cfg.permutation, 'yes'))
    error('multiple trials required in the data\n');
    
  elseif strcmp(cfg.permutation, 'yes')
    % compute the cross-spectral density matrix without resampling
    [dum, avg_aCf, avg_aCr, avg_aPr, avg_bCf, avg_bCr, avg_bPr] = prepare_resampled_data(cfg2 , aCf, aCr, aPr, bCf, bCr, bPr);
    % compute the cross-spectral density matrix with random permutation
    [dum, rnd_aCf, rnd_aCr, rnd_aPr, rnd_bCf, rnd_bCr, rnd_bPr] = prepare_resampled_data(cfg, aCf, aCr, aPr, bCf, bCr, bPr);
    % concatenate the different resamplings
    Cf = cat(1, reshape(avg_aCf, [1 Nchans Nchans]), reshape(avg_bCf, [1 Nchans Nchans]), rnd_aCf, rnd_bCf);
    Cr = cat(1, reshape(avg_aCr, [1 Nchans 1     ]), reshape(avg_bCr, [1 Nchans 1     ]), rnd_aCr, rnd_bCr);
    Pr = cat(1, reshape(avg_aPr, [1 1      1     ]), reshape(avg_bPr, [1 1      1     ]), rnd_aPr, rnd_bPr);
    % clear temporary working copies
    clear avg_aCf avg_aCr avg_aPr avg_bCf avg_bCr avg_bPr
    clear rnd_aCf rnd_aCr rnd_aPr rnd_bCf rnd_bCr rnd_bPr
    % the order of the resamplings should be [avgA avgB rndA rndB rndA rndB rndA rndB ....]
    Nrepetitions = 2*cfg.numpermutation + 2;
    order = [1 2 3:2:Nrepetitions 4:2:Nrepetitions];
    Cf = Cf(order,:,:);
    Cr = Cr(order,:,:);
    Pr = Pr(order,:,:);
    
  elseif strcmp(cfg.randomization, 'yes')
    % compute the cross-spectral density matrix without resampling
    [dum, avg_aCf, avg_aCr, avg_aPr, avg_bCf, avg_bCr, avg_bPr] = prepare_resampled_data(cfg2 , aCf, aCr, aPr, bCf, bCr, bPr);
    % compute the cross-spectral density matrix with random resampling
    [dum, rnd_aCf, rnd_aCr, rnd_aPr, rnd_bCf, rnd_bCr, rnd_bPr] = prepare_resampled_data(cfg, aCf, aCr, aPr, bCf, bCr, bPr);
    % concatenate the different resamplings
    Cf = cat(1, reshape(avg_aCf, [1 Nchans Nchans]), reshape(avg_bCf, [1 Nchans Nchans]), rnd_aCf, rnd_bCf);
    Cr = cat(1, reshape(avg_aCr, [1 Nchans 1     ]), reshape(avg_bCr, [1 Nchans 1     ]), rnd_aCr, rnd_bCr);
    Pr = cat(1, reshape(avg_aPr, [1 1      1     ]), reshape(avg_bPr, [1 1      1     ]), rnd_aPr, rnd_bPr);
    % clear temporary working copies
    clear avg_aCf avg_aCr avg_aPr avg_bCf avg_bCr avg_bPr
    clear rnd_aCf rnd_aCr rnd_aPr rnd_bCf rnd_bCr rnd_bPr
    % the order of the resamplings should be [avgA avgB rndA rndB rndA rndB rndA rndB ....]
    Nrepetitions = 2*cfg.numrandomization + 2;
    order = [1 2 3:2:Nrepetitions 4:2:Nrepetitions];
    Cf = Cf(order,:,:);
    Cr = Cr(order,:,:);
    Pr = Pr(order,:,:);
    
  elseif strcmp(cfg.jackknife, 'yes')
    % compute the cross-spectral density matrix with jackknife resampling
    [cfg, Cf, Cr, Pr] = prepare_resampled_data(cfg, Cf, Cr, Pr);
    Nrepetitions = Ntrials;
    
  elseif strcmp(cfg.bootstrap, 'yes')
    % compute the cross-spectral density matrix with bootstrap resampling
    [cfg, Cf, Cr, Pr] = prepare_resampled_data(cfg, Cf, Cr, Pr);
    Nrepetitions = cfg.numbootstrap;
    
  elseif strcmp(cfg.pseudovalue, 'yes')
    % compute the cross-spectral density matrix with pseudovalue resampling
    [cfg, Cf, Cr, Pr] = prepare_resampled_data(cfg, Cf, Cr, Pr);
    Nrepetitions = Ntrials+1;
    
  elseif strcmp(cfg.singletrial, 'yes')
    % The idea is that beamformer uses the average covariance to construct the
    % filter and applies it to the single trial covariance/csd. The problem
    % is that beamformer will use the averaged covariance/csd to estimate the
    % power and not the single trial covariance/csd
    error('this option contains a bug, and is therefore not supported at the moment');
    Cf = Cf; % FIXME, should be averaged and repeated for each trial
    Cr = Cr; % FIXME, should be averaged and repeated for each trial
    Pr = Pr; % FIXME, should be averaged and repeated for each trial
    Nrepetitions = Ntrials;
    
  elseif strcmp(cfg.rawtrial, 'yes')
    % keep all the individual trials, do not average them
    Cf = Cf;
    Cr = Cr;
    Pr = Pr;
    Nrepetitions = Ntrials;
    
  elseif Ntrials>1
    % compute the average from the individual trials
    Cf = reshape(sum(Cf, 1) / Ntrials, [Nchans Nchans]);
    Cr = reshape(sum(Cr, 1) / Ntrials, [Nchans 1]);
    Pr = reshape(sum(Pr, 1) / Ntrials, [1 1]);
    Nrepetitions = 1;
    
  elseif Ntrials==1
    % no rearrangement of trials is neccesary, the data already represents the average
    Cf = Cf;
    Cr = Cr;
    Pr = Pr;
    Nrepetitions = 1;
  end
  
  % reshape so that it also looks like one trial (out of many)
  if Nrepetitions==1
    Cf  = reshape(Cf , [1 Nchans Nchans]);
    Cr  = reshape(Cr , [1 Nchans 1]);
    Pr  = reshape(Pr , [1 1 1]);
  end
  
  % get the relevant low level options from the cfg and convert into key-value pairs
  tmpcfg = cfg.(cfg.method);
  % disable console feedback for the low-level function in case of multiple
  % repetitions
  if Nrepetitions > 1
    tmpcfg.feedback = 'none';
  end
  optarg = ft_cfg2keyval(tmpcfg);
  
  if Nrepetitions > 1
    ft_progress('init', cfg.(cfg.method).feedback, 'scanning repetition...');
  end
  
  for i=1:Nrepetitions
    size_Cf    = size(Cf);
    squeeze_Cf = reshape(Cf(i,:,:), size_Cf(2:end));    
    
    if Nrepetitions > 1
      ft_progress(i/Nrepetitions, 'scanning repetition %d from %d', i, Nrepetitions);
    end
    
    switch cfg.method
      case 'dics'
        if strcmp(submethod, 'dics_power')
          dip(i) = beamformer_dics(grid, sens, headmodel, [],  squeeze_Cf, optarg{:});
        elseif strcmp(submethod, 'dics_refchan')
          dip(i) = beamformer_dics(grid, sens, headmodel, [],  squeeze_Cf, optarg{:}, 'Cr', Cr(i,:), 'Pr', Pr(i));
        elseif strcmp(submethod, 'dics_refdip')
          dip(i) = beamformer_dics(grid, sens, headmodel, [],  squeeze_Cf, optarg{:}, 'refdip', cfg.refdip);
        end
      case 'pcc'
        if ~isempty(avg) && istrue(cfg.rawtrial)
          % FIXME added by jansch because an appropriate subselection of avg
          % should be done first (i.e. select the tapers that belong to this
          % repetition
          error('rawtrial in combination with pcc has been temporarily disabled');
        else
          dip(i) = beamformer_pcc(grid, sens, headmodel, avg, squeeze_Cf, optarg{:}, 'refdip', cfg.refdip, 'refchan', refchanindx, 'supdip', cfg.supdip, 'supchan', supchanindx);
        end
      case 'eloreta'
        dip(i) = ft_eloreta(grid, sens, headmodel, avg, squeeze_Cf, optarg{:});
      case 'mne'
        dip(i) = minimumnormestimate(grid, sens, headmodel, avg, optarg{:});
      case 'harmony'
        dip(i) = harmony(grid, sens, headmodel, avg, optarg{:});
        % error(sprintf('method ''%s'' is unsupported for source reconstruction in the frequency domain', cfg.method));
      case {'rv'}
        dip(i) = residualvariance(grid, sens, headmodel, avg, optarg{:}) ;
      case {'music'}
        error('method ''%s'' is currently unsupported for source reconstruction in the frequency domain', cfg.method);
      otherwise
    end
    
  end
  if Nrepetitions > 1
    ft_progress('close');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do time domain source reconstruction
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif istimelock && any(strcmp(cfg.method, {'lcmv', 'sam', 'mne','harmony', 'rv', 'music', 'pcc', 'mvl', 'sloreta', 'eloreta'}))
  
  % determine the size of the data
  Nsamples = size(data.avg,2);
  Nchans   = length(data.label);
  if isfield(data, 'cov') && length(size(data.cov))==3
    Ntrials = size(data.cov,1);
  elseif isfield(data, 'trial') && length(size(data.trial))==3
    Ntrials = size(data.trial,1);
  else
    Ntrials = 1;
  end
  
  if isfield(data, 'cov')
    % use the estimated data covariance matrix
    hascovariance = 1;
  else
    % add a identity covariance matrix, this simplifies the handling of the different source reconstruction methods
    % since the covariance is only used by some reconstruction methods and might not always be present in the data
    if Ntrials==1
      data.cov = eye(Nchans);
    else
      data.cov = zeros(Ntrials,Nchans,Nchans);
      for i=1:Ntrials
        data.cov(i,:,:) = eye(Nchans);
      end
    end
    hascovariance = 0;
    ft_warning('No covariance matrix found - will assume identity covariance matrix (mininum-norm solution)');
  end
  
  if strcmp(cfg.method, 'pcc')
    % HACK: requires some extra defaults
    if ~isfield(cfg, 'refdip'), cfg.refdip = []; end
    if ~isfield(cfg, 'supdip'), cfg.supdip = []; end
    
    % HACK: experimental code
    if hasbaseline
      error('not supported')
    end
    
    tmpcfg = [];
    tmpcfg.channel = cfg.channel(:)';
    if isfield(cfg, 'refchan')
      tmpcfg.channel = [tmpcfg.channel cfg.refchan(:)'];
    end
    if isfield(cfg, 'supchan')
      tmpcfg.channel = [tmpcfg.channel cfg.supchan(:)'];
    end
    
    % select the data in the channels of interest
    [dum, datchanindx] = match_str(tmpcfg.channel, data.label);
    if Ntrials==1
      data.avg = data.avg(datchanindx,:);
      data.cov = data.cov(datchanindx,datchanindx);
    else
      data.avg   = data.avg(datchanindx,:);
      data.cov   = data.cov(:,datchanindx,datchanindx);
      data.trial = data.trial(:,datchanindx,:);
    end
    data.label = data.label(datchanindx);
    
    if isfield(cfg, 'refchan') && ~isempty(cfg.refchan)
      [dum, refchanindx] = match_str(cfg.refchan, data.label);
    else
      refchanindx = [];
    end
    if isfield(cfg, 'supchan') && ~isempty(cfg.supchan)
      [dum, supchanindx] = match_str(cfg.supchan, data.label);
    else
      supchanindx = [];
    end
    Nchans = length(tmpcfg.channel); % update the number of channels
    
  else
    % HACK: use the default code
    % select the channels of interest
    [dum, datchanindx] = match_str(cfg.channel, data.label);
    if strcmp(data.dimord, 'chan_time')
      % It is in principle possible to have timelockanalysis with
      % keeptrial=yes and only a single trial in the raw data.
      % In that case the covariance should be represented as Nchan*Nchan
      data.avg = data.avg(datchanindx,:);
      %data.cov = reshape(data.cov, length(datchanindx), length(datchanindx));
      data.cov = data.cov(datchanindx,datchanindx);
    else
      data.avg   = data.avg(datchanindx,:);
      data.cov   = data.cov(:,datchanindx,datchanindx);
      data.trial = data.trial(:,datchanindx,:);
    end
    data.label = data.label(datchanindx);
    Nchans     = length(data.label);
  end
  
  if hasbaseline
    % baseline and active are only available together for resampling purposes,
    % hence I assume here that there are multiple trials in both
    baseline.avg   = baseline.avg(datchanindx,:);
    baseline.cov   = baseline.cov(:,datchanindx,datchanindx);
    baseline.trial = baseline.trial(:,datchanindx,:);
    % this is required for averaging 2 conditions using prepare_resampled_data
    cfg2 = [];
    cfg2.numcondition = 2;
  end
  
  % prepare the resampling of the trials, or average the data if multiple trials are present and no resampling is neccessary
  if (strcmp(cfg.jackknife, 'yes') || strcmp(cfg.bootstrap, 'yes') || strcmp(cfg.pseudovalue, 'yes') || strcmp(cfg.singletrial, 'yes') || strcmp(cfg.rawtrial, 'yes') || strcmp(cfg.randomization, 'yes')) && ~strcmp(data.dimord, 'rpt_chan_time')
    error('multiple trials required in the data\n');
    
  elseif strcmp(cfg.permutation, 'yes')
    % compute the average and covariance without resampling
    [dum, avgA, covA, avgB, covB] = prepare_resampled_data(cfg2 , data.trial, data.cov, baseline.trial, baseline.cov);
    % compute the average and covariance with random permutation
    [cfg, avRA, coRA, avRB, coRB] = prepare_resampled_data(cfg, data.trial, data.cov, baseline.trial, baseline.cov);
    % concatenate the different resamplings
    avg = cat(1, reshape(avgA, [1 Nchans Nsamples]), reshape(avgB, [1 Nchans Nsamples]), avRA, avRB);
    Cy  = cat(1, reshape(covA, [1 Nchans Nchans  ]), reshape(covB, [1 Nchans Nchans  ]), coRA, coRB);
    % clear temporary working copies
    clear avgA avgB covA covB
    clear avRA avRB coRA coRB
    % the order of the resamplings should be [avgA avgB randA randB randA randB randA randB ....]
    Nrepetitions = 2*cfg.numpermutation + 2;
    order = [1 2 3:2:Nrepetitions 4:2:Nrepetitions];
    avg = avg(order,:,:);
    Cy  = Cy (order,:,:);
    
  elseif strcmp(cfg.randomization, 'yes')
    % compute the average and covariance without resampling
    [dum, avgA, covA, avgB, covB] = prepare_resampled_data(cfg2 , data.trial, data.cov, baseline.trial, baseline.cov);
    % compute the average and covariance with random resampling
    [cfg, avRA, coRA, avRB, coRB] = prepare_resampled_data(cfg, data.trial, data.cov, baseline.trial, baseline.cov);
    % concatenate the different resamplings
    avg = cat(1, reshape(avgA, [1 Nchans Nsamples]), reshape(avgB, [1 Nchans Nsamples]), avRA, avRB);
    Cy  = cat(1, reshape(covA, [1 Nchans Nchans  ]), reshape(covB, [1 Nchans Nchans  ]), coRA, coRB);
    % clear temporary working copies
    clear avgA avgB covA covB
    clear avRA avRB coRA coRB
    % the order of the resamplings should be [avgA avgB randA randB randA randB randA randB ....]
    Nrepetitions = 2*cfg.numrandomization + 2;
    order = [1 2 3:2:Nrepetitions 4:2:Nrepetitions];
    avg = avg(order,:,:);
    Cy  = Cy (order,:,:);
    
  elseif strcmp(cfg.jackknife, 'yes')
    % compute the jackknife repetitions for the average and covariance
    [cfg, avg, Cy] = prepare_resampled_data(cfg, data.trial, data.cov);
    Nrepetitions = Ntrials;
    
  elseif strcmp(cfg.bootstrap, 'yes')
    % compute the bootstrap repetitions for the average and covariance
    [cfg, avg, Cy] = prepare_resampled_data(cfg, data.trial, data.cov);
    Nrepetitions = cfg.numbootstrap;
    
  elseif strcmp(cfg.pseudovalue, 'yes')
    % compute the pseudovalue repetitions for the average and covariance
    [cfg, avg, Cy] = prepare_resampled_data(cfg, data.trial, data.cov);
    Nrepetitions = Ntrials+1;
    
  elseif strcmp(cfg.singletrial, 'yes')
    % The idea is that beamformer uses the average covariance to construct the
    % filter and applies it to the single trial covariance/csd. The problem
    % is that beamformer will use the averaged covariance/csd to estimate the
    % power and not the single trial covariance/csd
    error('this option contains a bug, and is therefore not supported at the moment');
    % average the single-trial covariance matrices
    Cy = mean(data.cov,1);
    % copy the average covariance matrix for every individual trial
    Cy = repmat(Cy, [Ntrials 1 1]);
    % keep the single-trial ERFs, rename them to avg for convenience
    avg = data.trial;
    Nrepetitions = Ntrials;
    
  elseif strcmp(cfg.rawtrial, 'yes')
    % do not do any resampling, keep the single-trial covariances
    Cy = data.cov;
    % do not do any resampling, keep the single-trial ERFs (rename them to avg for convenience)
    avg = data.trial;
    Nrepetitions = Ntrials;
    
  elseif Ntrials>1
    % average the single-trial covariance matrices
    Cy  = reshape(mean(data.cov,1), [Nchans Nchans]);
    % select the average ERF
    avg = data.avg;
    Nrepetitions = 1;
    
  elseif Ntrials==1
    % select the average covariance matrix
    Cy  = data.cov;
    % select the average ERF
    avg = data.avg;
    Nrepetitions = 1;
  end
  
  % reshape so that it also looks like one trial (out of many)
  if Nrepetitions==1
    Cy  = reshape(Cy , [1 Nchans Nchans]);
    avg = reshape(avg, [1 Nchans Nsamples]);
  end
  
  % get the relevant low level options from the cfg and convert into key-value pairs
  optarg = ft_cfg2keyval(getfield(cfg, cfg.method));
  
  % This is the place to check for the consistency of the channel order in
  % the pre-computed leadfields/spatial filters, and to correct for it, if
  % necessary. This pertains to bugs 1746 and 3029.
  if isfield(grid, 'label') && (isfield(grid, 'leadfield') || isfield(grid, 'filter'))
    % match the channels in the leadfields/filters with those in the data
    [i1, i2] = match_str(cfg.channel, grid.label);
    if ~isequal(i2(:), (1:numel(grid.label))')
      if isfield(grid, 'leadfield')
        fprintf('\n\nSubselecting/reordering the channels in the precomputed leadfields\n\n');
        inside_indx = find(grid.inside);
        for k = inside_indx(:)'
          grid.leadfield{k} = grid.leadfield{k}(i2, :);
        end
      end
      if isfield(grid, 'filter')
        fprintf('\n\nSubselecting/reordering the channels in the precomputed filters\n\n');
        inside_indx = find(grid.inside);
        for k = inside_indx(:)'
          grid.filter{k} = grid.filter{k}(:, i2);
        end
      end
      grid.label = grid.label(i2);
    end
    if ~isequal(i1(:), (1:numel(cfg.channel))')
      % this is not so easy to deal with, throw an error
      error('There''s a mismatch between the number/order of channels in the data, with respect to the channels in the precomputed leadfield/filter. This is not easy to solve automatically. Please look into this.');
    end
  end
  
  size_avg = [size(avg) 1];
  size_Cy  = [size(Cy) 1];
  if strcmp(cfg.method, 'lcmv')% && ~isfield(grid, 'filter'),
    for i = 1:Nrepetitions
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      squeeze_Cy  = reshape(Cy(i,:,:), [size_Cy(2)  size_Cy(3)]);
      fprintf('scanning repetition %d\n', i);
      dip(i) = beamformer_lcmv(grid, sens, headmodel, squeeze_avg, squeeze_Cy, optarg{:});
    end
    
    % the following has been disabled since it turns out to be wrong (see
    % bugzilla bug 2395)
    %   elseif 0 && strcmp(cfg.method, 'lcmv')
    %     %don't loop over repetitions (slow), but reshape the input data to obtain single trial timecourses efficiently
    %     %in the presence of filters pre-computed on the average (or whatever)
    %     tmpdat = reshape(permute(avg,[2 3 1]),[size_avg(2) size_avg(3)*size_avg(1)]);
    %     tmpdip = beamformer_lcmv(grid, sens, headmodel, tmpdat, squeeze(mean(Cy,1)), optarg{:});
    %     tmpmom = tmpdip.mom{tmpdip.inside(1)};
    %     sizmom = size(tmpmom);
    %
    %     for i=1:length(tmpdip.inside)
    %       indx = tmpdip.inside(i);
    %       tmpdip.mom{indx} = permute(reshape(tmpdip.mom{indx}, [sizmom(1) size_avg(3) size_avg(1)]), [3 1 2]);
    %     end
    %     try, tmpdip = rmfield(tmpdip, 'pow'); end
    %     try, tmpdip = rmfield(tmpdip, 'cov'); end
    %     try, tmpdip = rmfield(tmpdip, 'noise'); end
    %     for i=1:Nrepetitions
    %       dip(i).pos     = tmpdip.pos;
    %       dip(i).inside  = tmpdip.inside;
    %       dip(i).outside = tmpdip.outside;
    %       dip(i).mom     = cell(1,size(tmpdip.pos,1));
    %       if isfield(tmpdip, 'ori')
    %         dip(i).ori   = cell(1,size(tmpdip.pos,1));
    %       end
    %       dip(i).cov     = cell(1,size(tmpdip.pos,1));
    %       dip(i).pow     = nan(size(tmpdip.pos,1),1);
    %       for ii=1:length(tmpdip.inside)
    %         indx             = tmpdip.inside(ii);
    %         tmpmom           = reshape(tmpdip.mom{indx}(i,:,:),[sizmom(1) size_avg(3)]);
    %         dip(i).mom{indx} = tmpmom;
    %         if isfield(tmpdip, 'ori')
    %           dip(i).ori{indx} = tmpdip.ori{indx};
    %         end
    %
    %         % the following recovers the single trial power and covariance, but
    %         % importantly the latency over which the power is defined is the
    %         % latency of the event-related field in the input and not the
    %         % latency of the covariance window, which can differ from the
    %         % former
    %         dip(i).cov{indx} = (tmpmom*tmpmom')./size_avg(3);
    %         if isempty(cfg.lcmv.powmethod) || strcmp(cfg.lcmv.powmethod, 'trace')
    %           dip(i).pow(indx) = trace(dip(i).cov{indx});
    %         else
    %           [tmpu,tmps,tmpv] = svd(dip(i).cov{indx});
    %           dip(i).pow(indx) = tmps(1);
    %         end
    %       end
    %     end
    
  elseif strcmp(cfg.method, 'sloreta')
    for i=1:Nrepetitions
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      squeeze_Cy  = reshape(Cy(i,:,:), [size_Cy(2)  size_Cy(3)]);
      fprintf('scanning repetition %d\n', i);
      dip(i) = ft_sloreta(grid, sens, headmodel, squeeze_avg, squeeze_Cy, optarg{:});
    end
    
  elseif strcmp(cfg.method, 'eloreta'),
    for i=1:Nrepetitions
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      squeeze_Cy  = reshape(Cy(i,:,:), [size_Cy(2)  size_Cy(3)]);
      fprintf('scanning repetition %d\n', i);
      dip(i) = ft_eloreta(grid, sens, headmodel, squeeze_avg, squeeze_Cy, optarg{:});
    end
  elseif strcmp(cfg.method, 'sam')
    for i=1:Nrepetitions
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      squeeze_Cy  = reshape(Cy(i,:,:), [size_Cy(2)  size_Cy(3)]);
      fprintf('scanning repetition %d\n', i);
      squeeze_avg=reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      dip(i) = beamformer_sam(grid, sens, headmodel, squeeze_avg, squeeze_Cy, optarg{:});
    end
  elseif strcmp(cfg.method, 'pcc')
    for i=1:Nrepetitions
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      squeeze_Cy  = reshape(Cy(i,:,:), [size_Cy(2)  size_Cy(3)]);
      fprintf('scanning repetition %d\n', i);
      squeeze_avg=reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      dip(i) = beamformer_pcc(grid, sens, headmodel, squeeze_avg, squeeze_Cy, optarg{:}, 'refdip', cfg.refdip, 'refchan', refchanindx, 'supchan', supchanindx);
    end
  elseif strcmp(cfg.method, 'mne')
    for i=1:Nrepetitions
      fprintf('estimating current density distribution for repetition %d\n', i);
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      if hascovariance
        squeeze_Cy  = reshape(Cy(i,:,:), [size_Cy(2)  size_Cy(3)]);
        dip(i) = minimumnormestimate(grid, sens, headmodel, squeeze_avg, optarg{:}, 'noisecov', squeeze_Cy);
      else
        dip(i) = minimumnormestimate(grid, sens, headmodel, squeeze_avg, optarg{:});
      end
    end
  elseif strcmp(cfg.method, 'harmony')
    for i=1:Nrepetitions
      fprintf('estimating current density distribution for repetition %d\n', i);
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      if hascovariance
        squeeze_Cy  = reshape(Cy(i,:,:), [size_Cy(2)  size_Cy(3)]);
        dip(i) = harmony(grid, sens, headmodel, squeeze_avg, optarg{:}, 'noisecov', squeeze_Cy);
      else
        dip(i) = harmony(grid, sens, headmodel, squeeze_avg, optarg{:});
      end
    end
  elseif strcmp(cfg.method, 'rv')
    for i=1:Nrepetitions
      fprintf('estimating residual variance at each grid point for repetition %d\n', i);
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      dip(i) = residualvariance(grid, sens, headmodel, squeeze_avg,      optarg{:});
    end
  elseif strcmp(cfg.method, 'music')
    for i=1:Nrepetitions
      fprintf('computing multiple signal classification for repetition %d\n', i);
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      if hascovariance
        squeeze_Cy  = reshape(Cy(i,:,:), [size_Cy(2)  size_Cy(3)]);
        dip(i) = music(grid, sens, headmodel, squeeze_avg, 'cov', squeeze_Cy, optarg{:});
      else
        dip(i) = music(grid, sens, headmodel, squeeze_avg,                    optarg{:});
      end
    end
  elseif strcmp(cfg.method, 'mvl')
    for i=1:Nrepetitions
      fprintf('estimating current density distribution for repetition %d\n', i);
      fns = fieldnames(cfg);
      optarg = cell(1,length(fns));
      n=1;
      for c=1:length(fns)
        optarg{n} = fns{c};
        optarg{n+1} = cfg.(fns{c});
        n=n+2;
      end
      squeeze_avg = reshape(avg(i,:,:),[size_avg(2) size_avg(3)]);
      dip(i) = mvlestimate(grid, sens, headmodel, squeeze_avg, optarg{:});
    end
  else
    error('method ''%s'' is unsupported for source reconstruction in the time domain', cfg.method);
  end
  
elseif iscomp
  error('the use of component data in ft_sourceanalysis is disabled for the time being: if you encounter this error message and you need this functionality please contact the FieldTrip development team');
else
  error('the specified method ''%s'' combined with the input data of type ''%s'' are not supported', cfg.method, ft_datatype(data));
end % if freq or timelock or comp data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up and collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(grid, 'dim')
  % the source reconstruction was perfomed on a regular 3d volume, remember the dimensions of the volume
  source.dim = grid.dim;
end

if isfield(grid, 'tri')
  % the source reconstruction was perfomed on a tesselated cortical sheet, remember the triangles
  source.tri = grid.tri;
end

if exist('grid', 'var')
  source = copyfields(grid, source, {'pos', 'inside', 'leadfield', 'leadfielddimord', 'label'});%, 'filter'});
end

if exist('dip', 'var')
  % the fields in the dip structure might be more recent than those in the grid structure
  source = copyfields(dip, source, {'pos', 'inside', 'leadfield', 'leadfielddimord', 'label'});%, 'filter'});
  
  % prevent duplication of these fields when copying the content of dip into source.avg or source.trial
  dip    = removefields(dip,       {'pos', 'inside', 'leadfield', 'leadfielddimord', 'label'});%, 'filter'});
  
  if istrue(cfg.(cfg.method).keepfilter) && isfield(dip(1), 'filter')
    for k = 1:numel(dip)
      dip(k).label        = sens.label;
      dip(k).filterdimord = '{pos}_ori_chan';
    end
  end
end

if ~istrue(cfg.keepleadfield)
  % remove the precomputed leadfields from the output source (if present)
  source = removefields(source, {'leadfield' 'leadfielddimord' 'label'});
end

% remove the precomputed leadfields from the cfg regardless of what keepleadfield is saying
% it should not be kept in cfg, since there it takes up too much space
cfg.grid = removefields(cfg.grid, {'leadfield' 'leadfielddimord' 'filter' 'filterdimord' 'label'});

if strcmp(cfg.jackknife, 'yes')
  source.method = 'jackknife';
  source.trial = dip;
  source.df = Ntrials;
elseif strcmp(cfg.bootstrap, 'yes')
  source.method = 'bootstrap';
  source.trial = dip;
  source.df = Ntrials;
elseif strcmp(cfg.pseudovalue, 'yes')
  source.method = 'pseudovalue';
  source.trial = dip;
elseif strcmp(cfg.singletrial, 'yes')
  source.method = 'singletrial';
  source.trial = dip;
  source.df = Ntrials;     % is this correct?
elseif strcmp(cfg.rawtrial, 'yes')
  source.method = 'rawtrial';
  source.trial = dip;
  source.df = Ntrials;     % is this correct?
elseif strcmp(cfg.randomization, 'yes')
  % assign the randomized resamplings to the output, keeping the special meaning of trial 1 and 2 in mind
  source.method = 'randomization';
  source.avgA = dip(1);
  source.avgB = dip(2);
  source.trialA = dip(1+2*(1:cfg.numrandomization));
  source.trialB = dip(2+2*(1:cfg.numrandomization));
elseif strcmp(cfg.permutation, 'yes')
  % assign the randomized resamplings to the output, keeping the special meaning of trial 1 and 2 in mind
  source.method = 'permutation';
  source.avgA = dip(1);
  source.avgB = dip(2);
  source.trialA = dip(1+2*(1:cfg.numpermutation));
  source.trialB = dip(2+2*(1:cfg.numpermutation));
elseif (strcmp(cfg.jackknife, 'yes') || strcmp(cfg.bootstrap, 'yes') || strcmp(cfg.pseudovalue, 'yes') || strcmp(cfg.singletrial, 'yes') || strcmp(cfg.rawtrial, 'yes')) && strcmp(cfg.keeptrials, 'yes')
  % keep the source reconstruction for each repeated or resampled trial
  source.trial = dip;
elseif exist('dip', 'var')
  % it looks like beamformer analysis was done on an average input, keep the average source reconstruction
  source.method = 'average';
  source.avg = dip;
else
  % apparently no computations were performed
end

% remember the trialinfo
if (strcmp(cfg.keeptrials, 'yes') || strcmp(cfg.method, 'pcc')) && isfield(data, 'trialinfo')
  source.trialinfo = data.trialinfo;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data baseline
ft_postamble provenance source
ft_postamble history    source
ft_postamble savevar    source
