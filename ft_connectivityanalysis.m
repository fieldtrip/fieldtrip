function [stat] = ft_connectivityanalysis(cfg, data)

% FT_CONNECTIVITYANALYSIS computes various measures of connectivity between
% MEG/EEG channels or between source-level signals.
%
% Use as
%   stat = ft_connectivityanalysis(cfg, data)
%   stat = ft_connectivityanalysis(cfg, timelock)
%   stat = ft_connectivityanalysis(cfg, freq)
%   stat = ft_connectivityanalysis(cfg, source)
% where the first input argument is a configuration structure (see below)
% and the second argument is the output of FT_PREPROCESSING,
% FT_TIMELOCKANLAYSIS, FT_FREQANALYSIS, FT_MVARANALYSIS or FT_SOURCEANALYSIS.
%
% The different connectivity metrics are supported only for specific
% datatypes (see below).
%
% The configuration structure has to contain
%   cfg.method  =  string, can be
%     'amplcorr',  amplitude correlation, support for freq and source data
%     'coh',       coherence, support for freq, freqmvar and source data.
%                  For partial coherence also specify cfg.partchannel, see below.
%                  For imaginary part of coherency or coherency also specify
%                  cfg.complex, see below.
%     'csd',       cross-spectral density matrix, can also calculate partial
%                  csds - if cfg.partchannel is specified, support for freq
%                  and freqmvar data
%     'dtf',       directed transfer function, support for freq and
%                  freqmvar data
%     'granger',   granger causality, support for freq and freqmvar data
%     'pdc',       partial directed coherence, support for freq and
%                  freqmvar data
%     'plv',       phase-locking value, support for freq and freqmvar data
%     'powcorr',   power correlation, support for freq and source data
%     'powcorr_ortho', power correlation with single trial
%                  orthogonalisation, support for source data
%     'ppc'        pairwise phase consistency
%     'psi',       phaseslope index, support for freq and freqmvar data
%     'wpli',      weighted phase lag index (signed one,
%                  still have to take absolute value to get indication of
%                  strength of interaction. Note: measure has positive
%                  bias. Use wpli_debiased to avoid this.
%     'wpli_debiased'  debiased weighted phase lag index
%                  (estimates squared wpli)
%     'wppc'       weighted pairwise phase consistency
%     'corr'       Pearson correlation, support for timelock or raw data
%
% Additional configuration options are
%   cfg.channel    = Nx1 cell-array containing a list of channels which are
%     used for the subsequent computations. This only has an effect when
%     the input data is univariate. See FT_CHANNELSELECTION
%   cfg.channelcmb = Nx2 cell-array containing the channel combinations on
%     which to compute the connectivity. This only has an effect when the
%     input data is univariate. See FT_CHANNELCOMBINATION
%   cfg.trials     = Nx1 vector specifying which trials to include for the
%     computation. This only has an effect when the input data contains
%     repetitions.
%   cfg.feedback   = string, specifying the feedback presented to the user.
%     Default is 'none'. See FT_PROGRESS
%
% For specific methods the cfg can also contain
%   cfg.partchannel = cell-array containing a list of channels that need to
%     be partialized out, support for method 'coh', 'csd', 'plv'
%   cfg.complex     = 'abs' (default), 'angle', 'complex', 'imag', 'real',
%     '-logabs', support for method 'coh', 'csd', 'plv'
%   cfg.removemean  = 'yes' (default), or 'no', support for method
%     'powcorr' and 'amplcorr'.
%   cfg.bandwidth   = scalar, (default = Rayleigh frequency), needed for
%			'psi', half-bandwidth of the integration across frequencies (in Hz)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_PREPROCESSING, FT_TIMELOCKANALYSIS, FT_FREQANALYSIS,
% FT_MVARANALYSIS, FT_SOURCEANALYSIS, FT_NETWORKANALYSIS.
%
% For the implemented methods, see also FT_CONNECTIVITY_CORR,
% FT_CONNECTIVITY_GRANGER, FT_CONNECTIVITY_PPC, FT_CONNECTIVITY_WPLI,
% FT_CONNECTIVITY_PDC, FT_CONNECTIVITY_DTF, FT_CONNECTIVITY_PSI

% Undocumented options:
%   cfg.refindx             =
%   cfg.jackknife           =
%   cfg.method              = 'mi';
%   cfg.granger.block       =
%   cfg.granger.conditional =
%
% Methods to be implemented
%                 'xcorr',     cross correlation function
%                 'di',        directionality index
%                 'spearman'   spearman's rank correlation

% Copyright (C) 2009, Jan-Mathijs Schoffelen, Andre Bastos, Martin Vinck, Robert Oostenveld
% Copyright (C) 2010-2011, Jan-Mathijs Schoffelen, Martin Vinck
% Copyright (C) 2012-2013, Jan-Mathijs Schoffelen
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

% FIXME it should be checked carefully whether the following works
% check if the input data is valid for this function
% data = ft_checkdata(data, 'datatype', {'raw', 'timelock', 'freq', 'source'});

% set the defaults
cfg.feedback    = ft_getopt(cfg, 'feedback', 'none');
cfg.channel     = ft_getopt(cfg, 'channel', 'all');
cfg.channelcmb  = ft_getopt(cfg, 'channelcmb', {'all' 'all'});
cfg.refindx     = ft_getopt(cfg, 'refindx', 'all');
cfg.trials      = ft_getopt(cfg, 'trials', 'all', 1);
cfg.complex     = ft_getopt(cfg, 'complex', 'abs');
cfg.jackknife   = ft_getopt(cfg, 'jackknife', 'no');
cfg.removemean  = ft_getopt(cfg, 'removemean', 'yes');
cfg.partchannel = ft_getopt(cfg, 'partchannel', '');
cfg.parameter   = ft_getopt(cfg, 'parameter', []);

hasjack = (isfield(data, 'method') && strcmp(data.method, 'jackknife')) || (isfield(data, 'dimord') && strcmp(data.dimord(1:6), 'rptjck'));
hasrpt  = (isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'rpt'))) || (isfield(data, 'avg') && isfield(data.avg, 'mom')) || (isfield(data, 'trial') && isfield(data.trial, 'mom')); % FIXME old-fashioned pcc data
dojack  = strcmp(cfg.jackknife, 'yes');
normrpt = 0; % default, has to be overruled e.g. in plv, because of single replicate normalisation
normpow = 1; % default, has to be overruled e.g. in csd,

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  data = ft_selectdata(tmpcfg, data);
  [cfg, data] = rollback_provenance(cfg, data);
end

% select channels/channelcombination of interest and set the cfg-options accordingly
if isfield(data, 'label'),
  selchan = cell(0, 1);
  if ~isempty(cfg.channelcmb) && ~isequal(cfg.channelcmb, {'all' 'all'}),
    tmpcmb = ft_channelcombination(cfg.channelcmb, data.label);
    tmpchan = unique(tmpcmb(:));
    cfg.channelcmb = ft_channelcombination(cfg.channelcmb, tmpchan, 1);
    selchan = [selchan;unique(cfg.channelcmb(:))];
  end
  
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  selchan = [selchan;cfg.channel];
  if ~isempty(cfg.partchannel)
    cfg.partchannel = ft_channelselection(cfg.partchannel, data.label);
    selchan = [selchan; cfg.partchannel];
  end
  tmpcfg = [];
  tmpcfg.channel = unique(selchan);
  data = ft_selectdata(tmpcfg, data);
elseif isfield(data, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
  if ~isempty(cfg.partchannel)
    error('partialisation is only possible without linearly indexed bivariate data');
  end
  if ~isempty(cfg.channelcmb),
    % FIXME do something extra here
  end
  % FIXME call selectdata
end

% FIXME check which methods require hasrpt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data bookkeeping - ensure that the input data is appropriate for the method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
needrpt = 1; % logical flag to specify whether (pseudo)-repetitions are required in the lower level connectivity function (can be singleton)
switch cfg.method
  case {'coh' 'csd'}
    if ~isempty(cfg.partchannel)
      if hasrpt && ~hasjack,
        warning('partialisation on single trial observations is not supported, removing trial dimension');
        try
          data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'}, 'cmbrepresentation', 'fullfast');
          inparam = 'crsspctrm';
          hasrpt = 0;
        catch
          error('partial coherence/csd is only supported for input allowing for a all-to-all csd representation');
        end
      else
        try
          data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'}, 'cmbrepresentation', 'full');
          inparam = 'crsspctrm';
        catch
          error('partial coherence/csd is only supported for input allowing for a all-to-all csd representation');
        end
      end
    else
      data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq' 'source'});
      inparam = 'crsspctrm';
    end
    
    if strcmp(cfg.method, 'csd'),
      normpow = 0;
      outparam = 'crsspctrm';
    elseif strcmp(cfg.method, 'coh'),
      outparam = 'cohspctrm';
    end
    
    dtype = ft_datatype(data);
    switch dtype
      case 'source'
        if isempty(cfg.refindx), error('indices of reference voxels need to be specified'); end
        % if numel(cfg.refindx)>1, error('more than one reference voxel is not yet supported'); end
      otherwise
    end
    % FIXME think of accommodating partial coherence for source data with only a few references
  case {'wpli'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'wplispctrm';
    debiaswpli = 0;
    if hasjack, error('to compute wpli, data should be in rpt format'); end
  case {'wpli_debiased'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'wpli_debiasedspctrm';
    debiaswpli = 1;
    if hasjack, error('to compute wpli, data should be in rpt format'); end
  case {'ppc'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'ppcspctrm';
    weightppc = 0;
    if hasjack, error('to compute ppc, data should be in rpt format'); end
  case {'wppc'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'wppcspctrm';
    weightppc = 1;
    if hasjack, error('to compute wppc, data should be in rpt format'); end
  case {'plv'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq' 'source'});
    inparam = 'crsspctrm';
    outparam = 'plvspctrm';
    normrpt = 1;
  case {'corr'}
    data = ft_checkdata(data, 'datatype', {'raw' 'timelock'});
    if isfield(data, 'cov')
      % it looks like a timelock with a cov, which is perfectly valid as input
      data = ft_checkdata(data, 'datatype', 'timelock');
    else
      % it does not have a cov, the covariance will be computed on the fly further down
      data = ft_checkdata(data, 'datatype', 'raw');
    end
    inparam = 'cov';
    outparam = cfg.method;
  case {'amplcorr' 'powcorr'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq' 'source'});
    dtype = ft_datatype(data);
    switch dtype
      case {'freq' 'freqmvar'}
        inparam = 'powcovspctrm';
      case 'source'
        inparam = 'powcov';
        if isempty(cfg.refindx), error('indices of reference voxels need to be specified'); end
        % if numel(cfg.refindx)>1, error('more than one reference voxel is not yet supported'); end
      otherwise
    end
    outparam = [cfg.method, 'spctrm'];
  case {'granger' 'instantaneous_causality' 'total_interdependence'}
    % create subcfg for the spectral factorization
    if ~isfield(cfg, 'granger')
      cfg.granger = [];
    end
    cfg.granger.conditional = ft_getopt(cfg.granger, 'conditional', 'no');
    cfg.granger.block       = ft_getopt(cfg.granger, 'block', []);
    if isfield(cfg, 'channelcmb'),
      cfg.granger.channelcmb = cfg.channelcmb;
      cfg = rmfield(cfg, 'channelcmb');
    end
    data = ft_checkdata(data, 'datatype', {'mvar' 'freqmvar' 'freq'});
    inparam = {'transfer', 'noisecov', 'crsspctrm'};
    if strcmp(cfg.method, 'granger'),                 outparam = 'grangerspctrm'; end
    if strcmp(cfg.method, 'instantaneous_causality'), outparam = 'instantspctrm'; end
    if strcmp(cfg.method, 'total_interdependence'),   outparam = 'totispctrm';    end
    
    % check whether the frequency bins are more or less equidistant
    dfreq = diff(data.freq)./mean(diff(data.freq));
    assert(all(dfreq>0.999) && all(dfreq<1.001), ['non equidistant frequency bins are not supported for method ',cfg.method]);
    
  case {'dtf' 'pdc'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'transfer';
    outparam = [cfg.method, 'spctrm'];
  case {'psi'}
    
    cfg.bandwidth = ft_getopt(cfg, 'bandwidth', []);
    cfg.normalize = ft_getopt(cfg, 'normalize', 'no');
    assert(~isempty(cfg.bandwidth), 'you need to supply cfg.bandwidth with ''psi'' as method');
    
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'psispctrm';
    
    % check whether the frequency bins are more or less equidistant
    dfreq = diff(data.freq)./mean(diff(data.freq));
    assert(all(dfreq>0.999) && all(dfreq<1.001), 'non equidistant frequency bins are not supported for method ''psi''');
    
  case {'powcorr_ortho'}
    data = ft_checkdata(data, 'datatype', {'source', 'freq'});
    % inparam = 'avg.mom';
    inparam = 'mom';
    outparam = 'powcorrspctrm';
  case {'mi'}
    % create the subcfg for the mutual information
    if ~isfield(cfg, 'mi'), cfg.mi = []; end
    cfg.mi.numbin = ft_getopt(cfg.mi, 'numbin', 10);
    cfg.mi.lags   = ft_getopt(cfg.mi, 'lags',   0);
    
    % what are the input requirements?
    data = ft_checkdata(data, 'datatype', {'raw' 'timelock' 'freq' 'source'});
    dtype = ft_datatype(data);
    if strcmp(dtype, 'timelock')
      if ~isfield(data, 'trial')
        inparam = 'avg';
      else
        inparam = 'trial';
      end
      hasrpt = (isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'rpt')));
    elseif strcmp(dtype, 'raw')
      inparam = 'trial';
      hasrpt  = 1;
    elseif strcmp(dtype, 'freq')
      inparam = 'something';
    else
      inparam = 'something else';
    end
    outparam = 'mi';
    needrpt  = 1;
  case {'di'}
    % wat eigenlijk?
  otherwise
    error('unknown method % s', cfg.method);
end

dtype = ft_datatype(data);

% FIXME throw an error if cfg.complex~='abs', and dojack==1
% FIXME throw an error if no replicates and cfg.method='plv'
% FIXME trial selection has to be implemented still

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data bookkeeping - check whether the required inparam is present in the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(~isfield(data, inparam)) || (isfield(data, 'crsspctrm') && (ischar(inparam) && strcmp(inparam, 'crsspctrm'))),
  if iscell(inparam)
    % in the case of multiple inparams, use the first one to check the
    % input data (e.g. checking for 'transfer' for requested granger)
    inparam = inparam{1};
  end
  
  switch dtype
    case {'freq' 'freqmvar'}
      if strcmp(inparam, 'crsspctrm')
        if isfield(data, 'fourierspctrm')
          [data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm', 'crsspctrm', dtype, 'cmb', cfg.channelcmb, 'keeprpt', normrpt);
        elseif strcmp(inparam, 'crsspctrm') && isfield(data, 'powspctrm')
          % if input data is old-fashioned, i.e. contains powandcsd
          [data, powindx, hasrpt] = univariate2bivariate(data, 'powandcsd', 'crsspctrm', dtype, 'cmb', cfg.channelcmb, 'keeprpt', normrpt);
        elseif isfield(data, 'labelcmb')
          powindx = labelcmb2indx(data.labelcmb);
        else
          powindx = [];
        end
      elseif strcmp(inparam, 'powcovspctrm')
        if isfield(data, 'powspctrm'),
          [data, powindx] = univariate2bivariate(data, 'powspctrm', 'powcovspctrm', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'cmb', cfg.channelcmb, 'sqrtflag', strcmp(cfg.method, 'amplcorr'));
        elseif isfield(data, 'fourierspctrm'),
          [data, powindx] = univariate2bivariate(data, 'fourierspctrm', 'powcovspctrm', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'cmb', cfg.channelcmb, 'sqrtflag', strcmp(cfg.method, 'amplcorr'));
        end
      elseif strcmp(inparam, 'transfer')
        if isfield(data, 'fourierspctrm')
          % FIXME this is fast but throws away the trial dimension, consider
          % a way to keep trial information if needed, but use the fast way
          % if possible
          data = ft_checkdata(data, 'cmbrepresentation', 'fullfast');
          hasrpt = 0;
        elseif isfield(data, 'powspctrm')
          data = ft_checkdata(data, 'cmbrepresentation', 'full');
        end
        
       % convert the inparam back to cell array in the case of granger
        if strcmp(cfg.method, 'granger') || strcmp(cfg.method, 'instantaneous_causality') || strcmp(cfg.method, 'total_interdependence')
          inparam = {'transfer' 'noisecov' 'crsspctrm'};
          tmpcfg  = ft_checkconfig(cfg, 'createsubcfg', {'granger'});
          optarg  = ft_cfg2keyval(tmpcfg.granger);
        else
          tmpcfg  = ft_checkconfig(cfg, 'createsubcfg', {cfg.method});
          optarg  = ft_cfg2keyval(tmpcfg.(cfg.method));
        end
        
        % compute the transfer matrix
        data   = ft_connectivity_csd2transfer(data, optarg{:});
      end
      
    case 'source'
      if ischar(cfg.refindx) && strcmp(cfg.refindx, 'all')
        cfg.refindx = 1:size(data.pos,1);
      elseif ischar(cfg.refindx)
        error('cfg.refindx should be a 1xN vector, or ''all''');
      end
      if strcmp(inparam, 'crsspctrm')
        [data, powindx, hasrpt] = univariate2bivariate(data, 'mom', 'crsspctrm', dtype, 'cmb', cfg.refindx, 'keeprpt', 0);
        % [data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm', 'crsspctrm', dtype, 0, cfg.refindx, [], 1);
      elseif strcmp(inparam, 'powcov')
        if isfield(data, 'pow')
          [data, powindx, hasrpt] = univariate2bivariate(data, 'pow', 'powcov', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'cmb', cfg.refindx, 'sqrtflag', strcmp(cfg.method, 'amplcorr'), 'keeprpt', 0);
        elseif isfield(data, 'mom')
          [data, powindx, hasrpt] = univariate2bivariate(data, 'mom', 'powcov', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'cmb', cfg.refindx, 'sqrtflag', strcmp(cfg.method, 'amplcorr'), 'keeprpt', 0);
        end
      end
      
    case 'comp'
      [data, powindx, hasrpt] = univariate2bivariate(data, 'trial', 'cov', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'cmb', cfg.channelcmb, 'sqrtflag', false, 'keeprpt', 1);
      
  end % switch dtype
  
elseif (isfield(data, 'crsspctrm') && (ischar(inparam) && strcmp(inparam, 'crsspctrm')))
  % this means that there is a sparse crsspctrm in the data
  
else
  powindx = [];
end % ensure that the bivariate measure exists

% do some additional work if single trial normalisation is required
% for example when plv needs to be computed
if normrpt && hasrpt,
  if strcmp(inparam, 'crsspctrm'),
    tmp = data.(inparam);
    nrpt = size(tmp, 1);
    ft_progress('init', cfg.feedback, 'normalising...');
    for k = 1:nrpt
      ft_progress(k/nrpt, 'normalising amplitude of replicate % d from % d to 1\n', k, nrpt);
      tmp(k, :, :, :, :) = tmp(k, :, :, :, :)./abs(tmp(k, :, :, :, :));
    end
    ft_progress('close');
    data.(inparam) = tmp;
  end
end

% convert the labels for the partialisation channels into indices
% do the same for the labels of the channels that should be kept
% convert the labels in the output to reflect the partialisation
if ~isempty(cfg.partchannel)
  allchannel = ft_channelselection(cfg.channel, data.label);
  pchanindx = match_str(allchannel, cfg.partchannel);
  kchanindx = setdiff(1:numel(allchannel), pchanindx);
  keepchn = allchannel(kchanindx);
  cfg.pchanindx = pchanindx;
  cfg.allchanindx = kchanindx;
  partstr = '';
  for k = 1:numel(cfg.partchannel)
    partstr = [partstr, '-', cfg.partchannel{k}];
  end
  for k = 1:numel(keepchn)
    keepchn{k} = [keepchn{k}, '\', partstr(2:end)];
  end
  data.label = keepchn; % update labels to remove the partialed channels
  % FIXME consider keeping track of which channels have been partialised
else
  cfg.pchanindx = [];
  cfg.allchanindx = [];
end

% check if jackknife is required
if hasrpt && dojack && hasjack,
  % do nothing
elseif hasrpt && dojack && ~(exist('debiaswpli', 'var') || exist('weightppc', 'var')),
  % compute leave-one-outs
  % assume the inparam(s) are well-behaved, i.e. they have the 'rpt'
  % dimension as the first dimension
  if iscell(inparam)
    for k = 1:numel(inparam)
      nrpt   = size(data.(inparam{k}),1);
      sumdat = sum(data.(inparam{k}),1);
      data.(inparam{k}) = (sumdat(ones(nrpt,1),:,:,:,:,:) - data.(inparam{k}))./(nrpt-1);
      clear sumdat;
    end
  else
    nrpt   = size(data.(inparam),1);
    sumdat = sum(data.(inparam),1);
    data.(inparam) = (sumdat(ones(nrpt,1),:,:,:,:,:) - data.(inparam))./(nrpt-1);
    clear sumdat;
  end
  hasjack = 1;
elseif hasrpt && ~(exist('debiaswpli', 'var') || exist('weightppc', 'var') || strcmp(cfg.method, 'powcorr_ortho'))% || needrpt)
  % create dof variable
  if isfield(data, 'dof')
    dof = data.dof;
  elseif isfield(data, 'cumtapcnt')
    dof = sum(data.cumtapcnt);
  end
  tmpcfg = [];
  tmpcfg.avgoverrpt = 'yes';
  data = ft_selectdata(tmpcfg, data);
  hasrpt = 0;
else
  % nothing required
end

% ensure that the first dimension is singleton if ~hasrpt
if ~hasrpt && needrpt
  if ischar(inparam)
    data.(inparam) = reshape(data.(inparam), [1 size(data.(inparam))]);
  else
    for k = 1:numel(inparam)
      data.(inparam{k}) = reshape(data.(inparam{k}), [1 size(data.(inparam{k}))]);
    end
  end
  if isfield(data, 'dimord')
    data.dimord = ['rpt_', data.dimord];
  elseif ~strcmp(dtype, 'raw')
    data.([inparam, 'dimord']) = ['rpt_', data.([inparam, 'dimord'])];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the desired connectivity metric by calling the appropriate ft_connectivity_XXX function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.method
  case 'coh'
    % coherence (unsquared), if cfg.complex = 'imag' imaginary part of coherency
    optarg = {'complex', cfg.complex, 'dimord', data.dimord, 'feedback', cfg.feedback, 'pownorm', normpow, 'hasjack', hasjack};
    if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'csd'
    % cross-spectral density (e.g. useful if partialisation is required)
    optarg = {'complex', cfg.complex, 'dimord', data.dimord, 'feedback', cfg.feedback, 'pownorm', normpow, 'hasjack', hasjack};
    if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case {'wpli' 'wpli_debiased'}
    % weighted pli or debiased weighted phase lag index.
    optarg = {'feedback', cfg.feedback, 'dojack', dojack, 'debias', debiaswpli};
    [datout, varout, nrpt] = ft_connectivity_wpli(data.(inparam), optarg{:});
    
  case {'wppc' 'ppc'}
    % weighted pairwise phase consistency or pairwise phase consistency
    optarg = {'feedback', cfg.feedback, 'dojack', dojack, 'weighted', weightppc};
    [datout, varout, nrpt] = ft_connectivity_ppc(data.(inparam), optarg{:});
    
  case 'plv'
    % phase locking value
    optarg = {'complex', cfg.complex, 'dimord', data.dimord, 'feedback', cfg.feedback, 'pownorm', normpow, 'hasjack', hasjack};
    if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'amplcorr'
    % amplitude correlation
    if isfield(data, 'dimord'),
      dimord = data.dimord;
    else
      dimord = data.([inparam, 'dimord']);
    end
    optarg = {'feedback', cfg.feedback, 'dimord', dimord, 'complex', 'real', 'pownorm', 1, 'pchanindx', [], 'hasjack', hasjack};
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'powcorr'
    % power correlation
    if isfield(data, 'dimord'),
      dimord = data.dimord;
    else
      dimord = data.([inparam, 'dimord']);
    end
    optarg = {'feedback', cfg.feedback, 'dimord', dimord, 'complex', 'real', 'pownorm', 1, 'pchanindx', [], 'hasjack', hasjack};
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case {'granger' 'instantaneous_causality' 'total_interdependence'}
    % granger causality
    if ft_datatype(data, 'freq') || ft_datatype(data, 'freqmvar'),
      if isfield(data, 'labelcmb') && ~istrue(cfg.granger.conditional),
        % multiple pairwise non-parametric transfer functions
        % linearly indexed
        
        % The following is very slow, one may make assumptions regarding
        % the order of the channels -> csd2transfer gives combinations in
        % quadruplets, where the first and fourth are auto-combinations,
        % and the second and third are cross-combinations
        % powindx = labelcmb2indx(data.labelcmb);
        %
        % The following is not needed anymore, because ft_connectivity_granger
        % relies on some hard-coded assumptions for the channel-pair ordering.
        % Otherwise it becomes just too slow.
        % powindx = zeros(size(data.labelcmb));
        % for k = 1:size(powindx, 1)/4
        % ix = ((k-1)*4+1):k*4;
        % powindx(ix, :) = [1 1;4 1;1 4;4 4] + (k-1)*4;
        % end
        
        powindx = [];
        
        if isfield(data, 'label'),
          % this field should be removed
          data = rmfield(data, 'label');
        end
        
      elseif isfield(data, 'labelcmb') && istrue(cfg.granger.conditional),
        % conditional (blockwise) needs linearly represented cross-spectra,
        % that have been produced by ft_connectivity_csd2transfer
        %
        % each row in Nx2 cell-array tmp refers to a comparison
        % tmp{k, 1} represents the ordered blocks
        % for the full trivariate model: the second element drives the
        % first element, while the rest is partialed out.
        % tmp{k, 2} represents the ordered blocks where the driving block
        % is left out
        
        
        blocks  = unique(data.blockindx);
        nblocks = numel(blocks);
        
        cnt = 0;
        for k = 1:nblocks
          for m = (k+1):nblocks
            cnt  = cnt+1;
            rest = setdiff(reshape(blocks,[1 numel(blocks)]), [k m]); % make sure to reshape blocks into 1xn vector
            tmp{cnt, 1} = [k m rest];
            tmp{cnt, 2} = [k   rest];
            newlabelcmb{cnt, 1} = data.block(m).name; % note the index swap: convention is driver in left column
            newlabelcmb{cnt, 2} = data.block(k).name;
            cnt  = cnt+1;
            tmp{cnt, 1} = [m k rest];
            tmp{cnt, 2} = [m   rest];
            newlabelcmb{cnt, 1} = data.block(k).name;
            newlabelcmb{cnt, 2} = data.block(m).name;
          end
        end
        [cmbindx, n] = blockindx2cmbindx(data.labelcmb, {data.label data.blockindx}, tmp);
        powindx.cmbindx = cmbindx;
        powindx.n = n;
        data.labelcmb = newlabelcmb;
        
        if isfield(data, 'label')
          % this field should be removed
          data = rmfield(data, 'label');
        end
        
      elseif isfield(cfg.granger, 'block') && ~isempty(cfg.granger.block)
        % blockwise granger
        for k = 1:numel(cfg.granger.block)
          %newlabel{k, 1} = cat(2, cfg.granger.block(k).label{:});
          newlabel{k,1}  = cfg.granger.block(k).name;
          powindx{k,1}   = match_str(data.label, cfg.granger.block(k).label);
        end
        data.label = newlabel;
      else
        powindx = [];
      end
      % fs = cfg.fsample; % FIXME do we really need this, or is this related to how noisecov is defined and normalised?
      if ~exist('powindx', 'var'), powindx = []; end
      if strcmp(cfg.method, 'granger'),                 methodstr = 'granger';      end
      if strcmp(cfg.method, 'instantaneous_causality'), methodstr = 'instantaneous'; end
      if strcmp(cfg.method, 'total_interdependence'),   methodstr = 'total';        end
      optarg = {'hasjack', hasjack, 'method', methodstr, 'powindx', powindx, 'dimord', data.dimord};
      [datout, varout, nrpt] = ft_connectivity_granger(data.transfer, data.noisecov, data.crsspctrm, optarg{:});
    else
      error('granger for time domain data is not yet implemented');
    end
    
  case 'dtf'
    % directed transfer function
    if isfield(data, 'labelcmb'),
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
    optarg = {'feedback', cfg.feedback, 'powindx', powindx, 'hasjack', hasjack};
    hasrpt = ~isempty(strfind(data.dimord, 'rpt'));
    if hasrpt,
      nrpt = size(data.(inparam), 1);
      datin = data.(inparam);
    else
      nrpt = 1;
      datin = reshape(data.(inparam), [1 size(data.(inparam))]);
    end
    [datout, varout, nrpt] = ft_connectivity_dtf(datin, optarg{:});
    
  case 'pdc'
    % partial directed coherence
    if isfield(data, 'labelcmb'),
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
    optarg = {'feedback', cfg.feedback, 'powindx', powindx, 'hasjack', hasjack};
    hasrpt = ~isempty(strfind(data.dimord, 'rpt'));
    if hasrpt,
      nrpt = size(data.(inparam), 1);
      datin = data.(inparam);
    else
      nrpt = 1;
      datin = reshape(data.(inparam), [1 size(data.(inparam))]);
    end
    [datout, varout, nrpt] = ft_connectivity_pdc(datin, optarg{:});
    
  case 'psi'
    % phase slope index
    nbin = nearest(data.freq, data.freq(1)+cfg.bandwidth)-1;
    
    optarg = {'feedback', cfg.feedback, 'dimord', data.dimord, 'nbin', nbin, 'normalize', cfg.normalize, 'hasrpt', hasrpt, 'hasjack', hasjack};
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_psi(data.(inparam), optarg{:});
    
  case 'powcorr_ortho'
    % Joerg Hipp's power correlation method
    optarg = {'refindx', cfg.refindx, 'tapvec', data.cumtapcnt};
    if isfield(data, 'mom')
      % this is expected to be a single frequency
      %dat    = cat(2, data.mom{data.inside}).';
      
      % HACK
      dimord = getdimord(data, 'mom');
      dimtok = tokenize(dimord, '_');
      posdim = find(strcmp(dimtok,'{pos}'));
      posdim = 4; % we concatenate across positions...
      rptdim = find(~cellfun('isempty',strfind(dimtok,'rpt')));
      rptdim = rptdim-1; % the posdim has to be taken into account...
      dat    = cat(4, data.mom{data.inside});
      dat    = permute(dat,[posdim rptdim setdiff(1:ndims(dat),[posdim rptdim])]);
      
      datout = ft_connectivity_powcorr_ortho(dat, optarg{:});
    elseif strcmp(data.dimord, 'rpttap_chan_freq')
      % loop over all frequencies
      [nrpttap, nchan, nfreq] = size(data.fourierspctrm);
      datout = cell(1, nfreq);
      for i=1:length(data.freq)
        dat       = reshape(data.fourierspctrm(:,:,i)', nrpttap, nchan).';
        datout{i} = ft_connectivity_powcorr_ortho(dat, optarg{:});
      end
      datout = cat(3, datout{:});
      % HACK otherwise I don't know how to inform the code further down about the dimord
      data.dimord = 'rpttap_chan_chan_freq';
    else
      error('unsupported data representation');
    end
    varout = [];
    nrpt = numel(data.cumtapcnt);
    
  case 'mi'
    % mutual information using the information breakdown toolbox
    % presence of the toolbox is checked in the low-level function
    
    if ~strcmp(dtype, 'raw') && (numel(cfg.mi.lags)>1 || cfg.mi.lags~=0),
      error('computation of lagged mutual information is only possible with ''raw'' data in the input');
    end
    
    switch dtype
      case 'raw'
        % ensure the lags to be in samples, not in seconds.
        cfg.mi.lags = round(cfg.mi.lags.*data.fsample);
        
        dat = catnan(data.trial, max(abs(cfg.mi.lags)));
        
        
        if ischar(cfg.refindx) && strcmp(cfg.refindx, 'all')
          outdimord = 'chan_chan';
        elseif numel(cfg.refindx)==1,
          outdimord = 'chan';
        else
          error('at present cfg.refindx should be either ''all'', or scalar');
        end
        if numel(cfg.mi.lags)>1
          data.time = cfg.mi.lags./data.fsample;
          outdimord = [outdimord,'_time'];
        else
          data = rmfield(data, 'time');
        end
        
      case 'timelock'
        dat = data.(inparam);
        dat = reshape(permute(dat, [2 3 1]), [size(dat, 2) size(dat, 1)*size(dat, 3)]);
        data = rmfield(data, 'time');
        if ischar(cfg.refindx) && strcmp(cfg.refindx, 'all')
          outdimord = 'chan_chan';
        elseif numel(cfg.refindx)==1,
          outdimord = 'chan';
        else
          error('at present cfg.refindx should be either ''all'', or scalar');
        end
        
        %data.dimord = 'chan_chan';
      case 'freq'
        error('not yet implemented');
      case 'source'
        % for the time being work with mom
        % dat = cat(2, data.mom{data.inside}).';
        dat = cat(1, data.mom{data.inside});
        % dat = abs(dat);
    end
    optarg = {'numbin', cfg.mi.numbin, 'lags', cfg.mi.lags, 'refindx', cfg.refindx};
    [datout] = ft_connectivity_mutualinformation(dat, optarg{:});
    varout = [];
    nrpt = [];
    
  case 'corr'
    % pearson's correlation coefficient
    optarg = {'dimord', getdimord(data, inparam), 'feedback', cfg.feedback, 'hasjack', hasjack};
    if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'xcorr'
    % cross-correlation function
    error('method %s is not yet implemented', cfg.method);
    
  case 'spearman'
    % spearman's rank correlation
    error('method %s is not yet implemented', cfg.method);
    
  case 'di'
    % directionality index
    error('method %s is not yet implemented', cfg.method);
    
  otherwise
    error('unknown method %s', cfg.method);
    
end % switch method

% remove the auto combinations if necessary -> FIXME this is granger specific and thus could move to ft_connectivity_granger
if (strcmp(cfg.method, 'granger') || strcmp(cfg.method, 'instantaneous_causality') || strcmp(cfg.method, 'total_interdependence')) && isfield(cfg, 'granger') && isfield(cfg.granger, 'sfmethod') && strcmp(cfg.granger.sfmethod, 'bivariate'),
  % remove the auto-combinations based on the order in the data
  switch dtype
    case {'freq' 'freqmvar'}
      keepchn = 1:size(datout, 1);
      keepchn = mod(keepchn, 4)==2 | mod(keepchn, 4)==3;
      datout = datout(keepchn, :, :, :, :);
      if ~isempty(varout),
        varout = varout(keepchn, :, :, :, :);
      end
      data.labelcmb = data.labelcmb(keepchn, :);
    case 'source'
      % not yet implemented
  end
end

if exist('powindx', 'var') && ~isempty(powindx),
  % based on powindx
  switch dtype
    case {'freq' 'freqmvar'}
      if isfield(data, 'labelcmb') && ~isstruct(powindx),
        keepchn = powindx(:, 1) ~= powindx(:, 2);
        datout = datout(keepchn, :, :, :, :);
        if ~isempty(varout),
          if all(size(varout)==size(nrpt))
            nrpt = nrpt(keepchn, :, :, :, :);
          end
          varout = varout(keepchn, :, :, :, :);
        end
        data.labelcmb = data.labelcmb(keepchn, :);
      end
    case 'source'
      nvox = size(unique(data.pos(:, 1:3), 'rows'), 1);
      ncmb = size(data.pos, 1)/nvox-1;
      remove = (powindx(:, 1) == powindx(:, 2)) & ((1:size(powindx, 1))' > nvox*ncmb);
      keepchn = ~remove;
      
      datout = datout(keepchn, :, :, :, :);
      if ~isempty(varout),
        varout = varout(keepchn, :, :, :, :);
      end
      inside = false(zeros(1, size(data.pos, 1)));
      inside(data.inside) = true;
      inside = inside(keepchn);
      data.inside = find(inside)';
      data.outside = find(inside==0)';
      data.pos = data.pos(keepchn, :);
  end % switch dtype
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the output structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch dtype
  case {'freq' 'freqmvar'},
    stat = [];
    if isfield(data, 'label'),
      stat.label = data.label;
    end
    if isfield(data, 'labelcmb'),
      stat.labelcmb = data.labelcmb;
      
      % ensure the correct dimord in case the input was 'powandcsd'
      data.dimord = strrep(data.dimord, 'chan_', 'chancmb_');
    end
    tok = tokenize(data.dimord, '_');
    dimord = '';
    for k = 1:numel(tok)
      if isempty(strfind(tok{k}, 'rpt'))
        dimord = [dimord, '_', tok{k}];
      end
    end
    stat.dimord = dimord(2:end);
    stat.(outparam) = datout;
    if ~isempty(varout),
      stat.([outparam, 'sem']) = (varout./nrpt).^0.5;
    end
    
  case 'timelock'
    stat = [];
    if isfield(data, 'label'),
      stat.label = data.label;
    end
    if isfield(data, 'labelcmb'),
      stat.labelcmb = data.labelcmb;
    end
    
    % deal with the dimord
    if exist('outdimord', 'var'),
      stat.dimord = outdimord;
    else
      % guess
      tok = tokenize(getdimord(data, inparam), '_');
      dimord = '';
      for k = 1:numel(tok)
        if isempty(strfind(tok{k}, 'rpt'))
          dimord = [dimord, '_', tok{k}];
        end
      end
      stat.dimord = dimord(2:end);
    end
    
    stat.(outparam) = datout;
    if ~isempty(varout),
      stat.([outparam, 'sem']) = (varout./nrpt).^0.5;
    end
    
  case 'source'
    stat = keepfields(data, {'pos', 'dim', 'transform', 'inside', 'outside'});
    stat.(outparam) = datout;
    if ~isempty(varout),
      stat.([outparam, 'sem']) = (varout/nrpt).^0.5;
    end
    
  case 'raw'
    stat = [];
    stat.label = data.label;
    stat.(outparam) = datout;
    if ~isempty(varout),
      stat.([outparam, 'sem']) = (varout/nrpt).^0.5;
    end
    if exist('outdimord', 'var'),
      stat.dimord = outdimord;
    end
end % switch dtype

if isfield(stat, 'dimord')
  dimtok = tokenize(stat.dimord, '_');
  % these dimensions in the output data must come from the input data
  if any(strcmp(dimtok, 'time')), stat.time = data.time; end
  if any(strcmp(dimtok, 'freq')), stat.freq = data.freq; end
else
  % just copy them over, alhtough we don't know for sure whether they are needed in the output
  if isfield(data, 'freq'), stat.freq = data.freq; end
  if isfield(data, 'time'), stat.time = data.time; end
end
if isfield(data, 'grad'), stat.grad = data.grad; end
if isfield(data, 'elec'), stat.elec = data.elec; end
if exist('dof', 'var'), stat.dof = dof; end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance stat
ft_postamble history    stat
ft_postamble savevar    stat


%-------------------------------------------------------------------------------
%subfunction to concatenate data with nans in between, needed for
%time-shifted mi
function [datamatrix] = catnan(datacells, nnans)

nchan = size(datacells{1}, 1);
nsmp  = cellfun('size',datacells,2);
nrpt  = numel(datacells);

%---initialize
datamatrix = nan(nchan, sum(nsmp) + nnans*(nrpt+1));

%---fill the matrix
for k = 1:nrpt
  if k==1,
    begsmp = 1+nnans;
    endsmp = nsmp(1)+nnans;
  else
    begsmp = k*nnans + sum(nsmp(1:(k-1))) + 1;
    endsmp = k*nnans + sum(nsmp(1:k));
  end
  datamatrix(:,begsmp:endsmp) = datacells{k};
end
