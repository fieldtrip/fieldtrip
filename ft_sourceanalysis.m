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
%                    'mvl'   multivariate Laplace source localization
% The DICS and PCC methods are for frequency domain data, all other methods
% are for time domain data.
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
% You can also use the FT_PREPARE_LEADFIELD function to create a grid with
% dipole positions and with precomputed leadfields. 
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
%   cfg.vol           = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%   cfg.hdmfile       = name of file containing the volume conduction model, see FT_READ_VOL
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.grad          = structure with gradiometer definition, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%   cfg.gradfile      = name of file containing the gradiometer definition, see FT_READ_SENS
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
% See also FT_SOURCEDESCRIPTIVES, FT_SOURCESTATISTICS, FT_PREPARE_LEADFIELD,
% FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

% Undocumented local options:
% cfg.numcomponents
% cfg.refchannel
% cfg.trialweight   = 'equal' or 'proportional'
% cfg.powmethod     = 'lambda1' or 'trace'
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
% cfg.mri
% cfg.smooth
% cfg.sourceunits
% cfg.threshold
% cfg.symmetry
%
% This function depends on FT_PREPARE_FREQ_MATRICES which has the following options:
% cfg.channel (default set in FT_SOURCEANALYSIS: cfg.channel = 'all'), documented
% cfg.dicsfix
% cfg.frequency, documented
% cfg.latency, documented
% cfg.refchan, documented
%
% This function depends on FT_PREPARE_RESAMPLED_DATA which has the following options:
% cfg.jackknife (default set in FT_SOURCEANALYSIS: cfg.jackknife = 'no'), documented
% cfg.numbootstrap, documented
% cfg.numcondition (set in FT_SOURCEANALYSIS: cfg.numcondition = 2)
% cfg.numpermutation (default set in FT_SOURCEANALYSIS: cfg.numpermutation = 100), documented
% cfg.numrandomization (default set in FT_SOURCEANALYSIS: cfg.numrandomization = 100), documented
% cfg.permutation (default set in FT_SOURCEANALYSIS: cfg.permutation = 'no'), documented
% cfg.pseudovalue (default set in FT_SOURCEANALYSIS: cfg.pseudovalue = 'no'), documented
% cfg.randomization (default set in FT_SOURCEANALYSIS: cfg.randomization = 'no'), documented
%
% This function depends on FT_PREPARE_VOL_SENS which has the following options:
% cfg.channel, (default set in FT_SOURCEANALYSIS: cfg.channel = 'all'), documented
% cfg.elec, documented
% cfg.elecfile, documented
% cfg.grad, documented
% cfg.gradfile, documented
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented
%
% This function depends on FT_PREPARE_LEADFIELD which has the following options:
% cfg.feedback, (default set in FT_SOURCEANALYSIS: cfg.feedback = 'text'), documented
% cfg.grid, documented
% cfg.lbex
% cfg.normalize (default set in FT_SOURCEANALYSIS), documented
% cfg.previous
% cfg.reducerank (default set in FT_SOURCEANALYSIS), documented
% cfg.sel50p
% cfg.version

% Copyright (c) 2003-2008, Robert Oostenveld, F.C. Donders Centre
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
ft_preamble init
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar data baseline

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'timelock', 'freq', 'comp'}, 'feedback', 'yes');
if nargin>2
  baseline = ft_checkdata(baseline, 'datatype', {'timelock', 'freq', 'comp'}, 'feedback', 'yes');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',     {'jacknife',   'jackknife'});
cfg = ft_checkconfig(cfg, 'renamed',     {'refchannel', 'refchan'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'power',           'dics'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'coh_refchan',     'dics'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'coh_refdip',      'dics'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'dics_cohrefchan', 'dics'});
cfg = ft_checkconfig(cfg, 'renamedval',  {'method', 'dics_cohrefdip',  'dics'});
cfg = ft_checkconfig(cfg, 'forbidden',   {'parallel', 'trials'});

% determine the type of input data
if isfield(data, 'freq')
  isfreq     = 1;
  iscomp     = 0;
  istimelock = 0;
elseif isfield(data, 'topo')
  isfreq     = 0;
  iscomp     = 1;
  istimelock = 0;
elseif isfield(data, 'time')
  iscomp     = 0;
  isfreq     = 0;
  istimelock = 1;
else
  error('input data is not recognized');
end

% set the defaults
if ~isfield(cfg, 'method') && istimelock, cfg.method = 'lcmv';      end
if ~isfield(cfg, 'method') && isfreq,     cfg.method = 'dics';      end
if ~isfield(cfg, cfg.method),             cfg.(cfg.method) = [];    end
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

% if ~isfield(cfg, 'reducerank'),     cfg.reducerank = 'no';      end  %
% the default for this depends on EEG/MEG and is set below
% put the low-level options pertaining to the source reconstruction method in their own field
% put the low-level options pertaining to the dipole grid in their own field

cfg = ft_checkconfig(cfg, 'createsubcfg',  {cfg.method, 'grid'});

cfg.(cfg.method).keepfilter    = ft_getopt(cfg.(cfg.method), 'keepfilter',    'no');
cfg.(cfg.method).keepcsd       = ft_getopt(cfg.(cfg.method), 'keepcsd',       'no');
cfg.(cfg.method).keepmom       = ft_getopt(cfg.(cfg.method), 'keepmom',       'yes');
cfg.(cfg.method).projectnoise  = ft_getopt(cfg.(cfg.method), 'projectnoise',  'no');
cfg.(cfg.method).feedback      = ft_getopt(cfg.(cfg.method), 'feedback',      'text');
cfg.(cfg.method).lambda        = ft_getopt(cfg.(cfg.method), 'lambda',        []);
cfg.(cfg.method).powmethod     = ft_getopt(cfg.(cfg.method), 'powmethod',     []);
cfg.(cfg.method).normalize     = ft_getopt(cfg.(cfg.method), 'normalize',     'no');

convertfreq = 0;
convertcomp = 0;
if ~istimelock && (strcmp(cfg.method, 'mne') || strcmp(cfg.method, 'rv') || strcmp(cfg.method, 'music'))
  % these timelock methods are also supported for frequency or component data
  if isfreq
    [data, cfg] = freq2timelock(cfg, data);
    convertfreq = 1;  % flag indicating that the data was converted
  elseif iscomp
    [data, cfg] = comp2timelock(cfg, data);
    convertcomp = 1;  % flag indicating that the data was converted
  end
  istimelock = 1;     % from now on the data can be treated as timelocked
  isfreq     = 0;
  iscomp     = 0;
end

% select only those channels that are present in the data
cfg.channel = ft_channelselection(cfg.channel, data.label);

if nargin>2 && (strcmp(cfg.randomization, 'no') && strcmp(cfg.permutation, 'no'))
  error('input of two conditions only makes sense if you want to randomize or permute');
elseif nargin<3 && (strcmp(cfg.randomization, 'yes') || strcmp(cfg.permutation, 'yes'))
  error('randomization or permutation requires that you give two conditions as input');
end

if isfield(cfg, 'latency') && istimelock
  error('specification of cfg.latency is only required for time-frequency data');
end

if sum([strcmp(cfg.jackknife, 'yes'), strcmp(cfg.bootstrap, 'yes'), strcmp(cfg.pseudovalue, 'yes'), strcmp(cfg.singletrial, 'yes'), strcmp(cfg.rawtrial, 'yes'), strcmp(cfg.randomization, 'yes'), strcmp(cfg.permutation, 'yes')])>1
  error('jackknife, bootstrap, pseudovalue, singletrial, rawtrial, randomization and permutation are mutually exclusive');
end

if strcmp(cfg.rawtrial,'yes') && isfield(cfg,'grid') && ~isfield(cfg.grid,'filter')
  error('Using each trial to compute its own filter is not currently recommended. Use this option only with precomputed filters in grid.filter');
end

if isfreq
  if  ~strcmp(data.dimord, 'chan_freq')          && ...
      ~strcmp(data.dimord, 'chan_freq_time')     && ...
      ~strcmp(data.dimord, 'rpt_chan_freq')      && ...
      ~strcmp(data.dimord, 'rpt_chan_freq_time') && ...
      ~strcmp(data.dimord, 'rpttap_chan_freq')   && ...
      ~strcmp(data.dimord, 'chancmb_freq')       && ...
      ~strcmp(data.dimord, 'rpt_chancmb_freq')   && ...
      ~strcmp(data.dimord, 'rpttap_chancmb_freq')  && ...
      ~strcmp(data.dimord, 'rpttap_chan_freq_time')
    error('dimord of input frequency data is not recognized');
  end
end

% collect and preprocess the electrodes/gradiometer and head model
[vol, sens, cfg] = prepare_headmodel(cfg, data);

% It might be that the number of channels in the data, the number of
% channels in the electrode/gradiometer definition and the number of
% channels in the localspheres volume conduction model are different.
% Hence a subset of the data channels will be used.
Nchans = length(cfg.channel);

% set the default for reducing the rank of the leadfields, this is an
% option to the specific method and will be passed on to the low-level
% function
if ~isfield(cfg.(cfg.method), 'reducerank')
  if ft_senstype(sens, 'meg')
    cfg.(cfg.method).reducerank = 2;
  else
    cfg.(cfg.method).reducerank = 3;
  end
end

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
end

if isfield(cfg.grid, 'filter')
  if numel(cfg.grid.filter) == size(grid.pos, 1)
    grid.filter = cfg.grid.filter;
  else
    warning_once('ignoring predefined filter as it does not match the grid''s dimension');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do frequency domain source reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfreq && any(strcmp(cfg.method, {'dics', 'pcc'}))
  
  if strcmp(cfg.method, 'pcc')
    % HACK: requires some extra defaults
    if ~isfield(cfg, 'refdip'), cfg.refdip = []; end
    if ~isfield(cfg, 'supdip'), cfg.supdip = []; end
    if ~isfield(cfg, 'refchan'), cfg.refchan = []; end
    if ~isfield(cfg, 'supchan'), cfg.supchan = []; end
    cfg.refchan = ft_channelselection(cfg.refchan, data.label);
    cfg.supchan = ft_channelselection(cfg.supchan, data.label);
    
    % HACK: use some experimental code
    if nargin>2,
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
    
    % if the input data has a complete fourier spectrum, project it through the filters
    % FIXME it was incorrect , since the
    % ' leads to a conjugate transposition check this in beamformer_pcc
    if isfield(data, 'fourierspctrm')
      [dum, datchanindx] = match_str(tmpcfg.channel, data.label);
      fbin = nearest(data.freq, cfg.frequency);
      if strcmp(data.dimord, 'chan_freq')
        avg = data.fourierspctrm(datchanindx, fbin);
      elseif strcmp(data.dimord, 'rpt_chan_freq') || strcmp(data.dimord, 'rpttap_chan_freq'),
        %avg = data.fourierspctrm(:, datchanindx, fbin)';
        avg = transpose(data.fourierspctrm(:, datchanindx, fbin));
      elseif strcmp(data.dimord, 'chan_freq_time')
        tbin = nearest(data.time, cfg.latency);
        avg = data.fourierspctrm(datchanindx, fbin, tbin);
      elseif strcmp(data.dimord, 'rpt_chan_freq_time') || strcmp(data.dimord, 'rpttap_chan_freq_time'),
        tbin = nearest(data.time, cfg.latency);
        %avg = data.fourierspctrm(:, datchanindx, fbin, tbin)';
        avg = transpose(data.fourierspctrm(:, datchanindx, fbin, tbin));
      end
    else
      avg = [];
    end
    
  else
    % HACK: use the default code
    % convert the input data, so that Cf, Cr and Pr contain either the average over all trials (if Ntrials==1)
    % or the individual cross-spectral-densities/powers of each individual trial (if Ntrials>1)
    [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices(cfg, data);
  end
  
  if strcmp(cfg.method, 'dics')
    % assign a descriptive name to each of the dics sub-methods, the default is power only
    if strcmp(cfg.method, 'dics') && isfield(cfg, 'refdip') && ~isempty(cfg.refdip);
      submethod = 'dics_refdip';
    elseif strcmp(cfg.method, 'dics') && isfield(cfg, 'refchan') && ~isempty(cfg.refchan);
      submethod = 'dics_refchan';
    else
      submethod = 'dics_power';
    end
  end
  
  % fill these with NaNs, so that I dont have to treat them separately
  if isempty(Cr), Cr = nan(Ntrials, Nchans, 1); end
  if isempty(Pr), Pr = nan(Ntrials, 1, 1); end
  
  if nargin>2
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
    ft_progress('init', 'text');
  end
  for i=1:Nrepetitions
    
    if Nrepetitions > 1
      ft_progress(i/Nrepetitions, 'scanning repetition %d of %d', i, Nrepetitions);
    end
    
    if     strcmp(cfg.method, 'dics') && strcmp(submethod, 'dics_power')
      dip(i) = beamformer_dics(grid, sens, vol, [],  squeeze(Cf(i,:,:)), optarg{:});
    elseif strcmp(cfg.method, 'dics') && strcmp(submethod, 'dics_refchan')
      dip(i) = beamformer_dics(grid, sens, vol, [],  squeeze(Cf(i,:,:)), optarg{:}, 'Cr', Cr(i,:), 'Pr', Pr(i));
    elseif strcmp(cfg.method, 'dics') && strcmp(submethod, 'dics_refdip')
      dip(i) = beamformer_dics(grid, sens, vol, [],  squeeze(Cf(i,:,:)), optarg{:}, 'refdip', cfg.refdip);
    elseif strcmp(cfg.method, 'pcc')
      dip(i) = beamformer_pcc(grid, sens, vol, avg, squeeze(Cf(i,:,:)), optarg{:}, 'refdip', cfg.refdip, 'refchan', refchanindx, 'supdip', cfg.supdip, 'supchan', supchanindx);
    else
      error(sprintf('method ''%s'' is unsupported for source reconstruction in the frequency domain', cfg.method));
    end
  end
  if Nrepetitions > 1
    ft_progress('close');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do time domain source reconstruction
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif istimelock && any(strcmp(cfg.method, {'lcmv', 'sam', 'mne', 'rv', 'music', 'pcc', 'mvl'}))
  
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
    warning_once('No covariance matrix found - will assume identity covariance matrix (mininum-norm solution)');
  end
  
  if strcmp(cfg.method, 'pcc')
    % HACK: requires some extra defaults
    if ~isfield(cfg, 'refdip'), cfg.refdip = []; end
    if ~isfield(cfg, 'supdip'), cfg.supdip = []; end
    
    % HACK: experimental code
    if nargin>2
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
  
  if nargin>2
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
  
  siz=[size(avg) 1];
  if strcmp(cfg.method, 'lcmv') && ~isfield(grid, 'filter'),
    for i=1:Nrepetitions
      squeeze_avg=reshape(avg(i,:,:),[siz(2) siz(3)]);
      fprintf('scanning repetition %d\n', i);
      dip(i) = beamformer_lcmv(grid, sens, vol, squeeze_avg, squeeze(Cy(i,:,:)), optarg{:});
    end
  elseif strcmp(cfg.method, 'lcmv')
    %don't loop over repetitions (slow), but reshape the input data to obtain single trial timecourses efficiently
    %in the presence of filters pre-computed on the average (or whatever)
    tmpdat = reshape(permute(avg,[2 3 1]),[siz(2) siz(3)*siz(1)]);
    tmpdip = beamformer_lcmv(grid, sens, vol, tmpdat, squeeze(mean(Cy,1)), optarg{:});
    tmpmom = tmpdip.mom{tmpdip.inside(1)};
    sizmom = size(tmpmom);

    for i=1:length(tmpdip.inside)
      indx = tmpdip.inside(i);
      tmpdip.mom{indx} = permute(reshape(tmpdip.mom{indx}, [sizmom(1) siz(3) siz(1)]), [3 1 2]);
    end
    try, tmpdip = rmfield(tmpdip, 'pow'); end
    try, tmpdip = rmfield(tmpdip, 'cov'); end
    try, tmpdip = rmfield(tmpdip, 'noise'); end
    for i=1:Nrepetitions
      dip(i).pos     = tmpdip.pos;
      dip(i).inside  = tmpdip.inside;
      dip(i).outside = tmpdip.outside;
      dip(i).mom     = cell(1,size(tmpdip.pos,1));
      if isfield(tmpdip, 'ori')
        dip(i).ori   = cell(1,size(tmpdip.pos,1));
      end
      dip(i).cov     = cell(1,size(tmpdip.pos,1));
      dip(i).pow     = nan(size(tmpdip.pos,1),1);
      for ii=1:length(tmpdip.inside)
        indx             = tmpdip.inside(ii);
        tmpmom           = reshape(tmpdip.mom{indx}(i,:,:),[sizmom(1) siz(3)]);
        dip(i).mom{indx} = tmpmom;
        if isfield(tmpdip, 'ori')
          dip(i).ori{indx} = tmpdip.ori{indx};
        end
       
        % the following recovers the single trial power and covariance, but
        % importantly the latency over which the power is defined is the
        % latency of the event-related field in the input and not the
        % latency of the covariance window, which can differ from the
        % former
        dip(i).cov{indx} = (tmpmom*tmpmom')./siz(3);
        if isempty(cfg.lcmv.powmethod) || strcmp(cfg.lcmv.powmethod, 'trace')
          dip(i).pow(indx) = trace(dip(i).cov{indx});
        else
          [tmpu,tmps,tmpv] = svd(dip(i).cov{indx});
          dip(i).pow(indx) = tmps(1);
        end
      end
    end
  elseif strcmp(cfg.method, 'sam')
    for i=1:Nrepetitions
      fprintf('scanning repetition %d\n', i);
      squeeze_avg=reshape(avg(i,:,:),[siz(2) siz(3)]);
      dip(i) = beamformer_sam(grid, sens, vol, squeeze_avg, squeeze(Cy(i,:,:)), optarg{:});
    end
  elseif strcmp(cfg.method, 'pcc')
    for i=1:Nrepetitions
      fprintf('scanning repetition %d\n', i);
      squeeze_avg=reshape(avg(i,:,:),[siz(2) siz(3)]);
      dip(i) = beamformer_pcc(grid, sens, vol, squeeze_avg, squeeze(Cy(i,:,:)), optarg{:}, 'refdip', cfg.refdip, 'refchan', refchanindx, 'supchan', supchanindx);
    end
  elseif strcmp(cfg.method, 'mne')
    for i=1:Nrepetitions
      fprintf('estimating current density distribution for repetition %d\n', i);
      squeeze_avg=reshape(avg(i,:,:),[siz(2) siz(3)]);
      if hascovariance
        dip(i) = minimumnormestimate(grid, sens, vol, squeeze_avg, optarg{:}, 'noisecov', squeeze(Cy(i,:,:)));
      else
        dip(i) = minimumnormestimate(grid, sens, vol, squeeze_avg, optarg{:});
      end
    end
  elseif strcmp(cfg.method, 'rv')
    for i=1:Nrepetitions
      fprintf('estimating residual variance at each grid point for repetition %d\n', i);
      squeeze_avg=reshape(avg(i,:,:),[siz(2) siz(3)]);
      dip(i) = residualvariance(grid, sens, vol, squeeze_avg,      optarg{:});
    end
  elseif strcmp(cfg.method, 'music')
    for i=1:Nrepetitions
      fprintf('computing multiple signal classification for repetition %d\n', i);
      squeeze_avg=reshape(avg(i,:,:),[siz(2) siz(3)]);
      if hascovariance
        dip(i) = music(grid, sens, vol, squeeze_avg, 'cov', squeeze(Cy(i,:,:)), optarg{:});
      else
        dip(i) = music(grid, sens, vol, squeeze_avg,                            optarg{:});
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
      squeeze_avg=reshape(avg(i,:,:),[siz(2) siz(3)]);
      dip(i) = mvlestimate(grid, sens, vol, squeeze_avg, optarg{:});
    end
  else
    error(sprintf('method ''%s'' is unsupported for source reconstruction in the time domain', cfg.method));
  end
  
end % if freq or timelock data

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

if istimelock
  % add the time axis to the output
  source.time = data.time;
elseif iscomp
  % FIXME, add the component numbers to the output
elseif isfreq
  % add the frequency axis to the output
  cfg.frequency    = data.freq(nearest(data.freq, cfg.frequency));
  source.freq = cfg.frequency;
  if isfield(data, 'time') && isfield(cfg, 'latency')
    cfg.latency    = data.time(nearest(data.time, cfg.latency));
    source.time    = cfg.latency;
  end
  if isfield(data, 'cumtapcnt'),
    source.cumtapcnt = data.cumtapcnt;
  end
end

if exist('dip', 'var')
  % do some cleaning up, keep the dipole positions etc. in the global structure and not in each trial
  source.pos     = dip(1).pos;
  source.inside  = dip(1).inside;
  source.outside = dip(1).outside;
  dip = rmfield(dip, 'pos');
  dip = rmfield(dip, 'inside');
  dip = rmfield(dip, 'outside');
  if isfield(dip(1), 'leadfield')
    source.leadfield = dip(1).leadfield;
    dip = rmfield(dip, 'leadfield');
  end
elseif exist('grid', 'var')
  % no scanning has been done, probably only the leadfield has been computed
  try, source.pos       = grid.pos;       end
  try, source.inside    = grid.inside;    end
  try, source.outside   = grid.outside;   end
  try, source.leadfield = grid.leadfield; end
end

if strcmp(cfg.keepleadfield, 'yes') && ~isfield(source, 'leadfield')
  % add the precomputed leadfields to the output source
  source.leadfield = grid.leadfield;
elseif strcmp(cfg.keepleadfield, 'no') && isfield(source, 'leadfield')
  % remove the precomputed leadfields from the output source
  source = rmfield(source, 'leadfield');
end


% remove the precomputed leadfields from the cfg regardless of what keepleadfield is saying
% it should not be kept in cfg, since there it takes up too much space
if isfield(cfg, 'grid') && isfield(cfg.grid, 'leadfield')
  cfg.grid = rmfield(cfg.grid, 'leadfield');
end

if convertfreq
  % FIXME, convert the source reconstruction back to a frequency representation
elseif convertcomp
  % FIXME, convert the source reconstruction back to a component representation
end

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
elseif exist('dip', 'var')
  % it looks like beamformer analysis was done on an average input, keep the average source reconstruction
  source.method = 'average';
  source.avg = dip;
else
  % apparently no computations were performed
end

if (strcmp(cfg.jackknife, 'yes') || strcmp(cfg.bootstrap, 'yes') || strcmp(cfg.pseudovalue, 'yes') || strcmp(cfg.singletrial, 'yes') || strcmp(cfg.rawtrial, 'yes')) && strcmp(cfg.keeptrials, 'yes')
  % keep the source reconstruction for each repeated or resampled trial
  source.trial  = dip;
end

% remember the trialinfo
if (strcmp(cfg.keeptrials, 'yes') || strcmp(cfg.method, 'pcc')) && isfield(data, 'trialinfo')
  source.trialinfo = data.trialinfo;
end

% clean up the output source structure, e.g. ensure that the size of output
% elements is npos*1 and not 1*npos
source = ft_checkdata(source);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
if nargin==2
  ft_postamble previous data
elseif nargin==3
  ft_postamble previous data baseline
end
ft_postamble history source
ft_postamble savevar source
