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
%     'dtf',       directed transfer function, support for freq and freqmvar data
%     'granger',   granger causality, support for freq and freqmvar data
%     'pdc',       partial directed coherence, support for freq and freqmvar data
%     'plv',       phase-locking value, support for freq and freqmvar data
%     'powcorr',   power correlation, support for freq and source data
%     'powcorr_ortho', power correlation with single trial
%                  orthogonalisation, support for source data
%     'ppc'        pairwise phase consistency
%     'psi',       phaseslope index, support for freq and freqmvar data
%     'wpli',      weighted phase lag index (signed one, still have to
%                  take absolute value to get indication of strength of
%                  interaction. Note that this measure has a positive
%                  bias. Use wpli_debiased to avoid this.
%     'wpli_debiased'  debiased weighted phase lag index (estimates squared wpli)
%     'wppc'       weighted pairwise phase consistency
%     'corr'       Pearson correlation, support for timelock or raw data
%     'laggedcoherence', lagged coherence estimate
%
% Additional configuration options are
%   cfg.channel    = Nx1 cell-array containing a list of channels which are
%                    used for the subsequent computations. This only has an effect
%                    when the input data is univariate. See FT_CHANNELSELECTION

%   cfg.channelcmb = Nx2 cell-array containing the channel combinations on
%                    which to compute the connectivity. This only has an effect when
%                    the input data is univariate. See FT_CHANNELCOMBINATION
%   cfg.trials     = Nx1 vector specifying which trials to include for the
%                    computation. This only has an effect when the input data
%                    contains repetitions.
%   cfg.feedback   = string, specifying the feedback presented to the user. Default
%                    is 'none'. See FT_PROGRESS
%
% For specific methods the configuration can also contain
%   cfg.partchannel = cell-array containing a list of channels that need to be
%                     partialized out, support for method 'coh', 'csd', 'plv'
%   cfg.complex     = string, 'abs' (default), 'angle', 'complex', 'imag', 'real',
%                     '-logabs', support for method 'coh', 'csd', 'plv'
%   cfg.removemean  = string, 'yes' (default), or 'no', support for
%                     method 'powcorr' and 'amplcorr'.
%   cfg.bandwidth   = scalar, needed for 'psi', half-bandwidth of the integration
%                     across frequencies (in Hz, default is the Rayleigh frequency)
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
%   cfg.method              = 'mi'/'di'/'dfi';
%   cfg.granger.block       =
%   cfg.granger.conditional =
%
% Methods to be implemented
%                 'xcorr',     cross correlation function
%                 'di',        directionality index
%                 'spearman'   spearman's rank correlation

% Copyright (C) 2009, Jan-Mathijs Schoffelen, Andre Bastos, Martin Vinck, Robert Oostenveld
% Copyright (C) 2010-2011, Jan-Mathijs Schoffelen, Martin Vinck
% Copyright (C) 2012-2019, Jan-Mathijs Schoffelen
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
cfg.feedback    = ft_getopt(cfg, 'feedback',   'none');
cfg.channel     = ft_getopt(cfg, 'channel',    'all');
cfg.channelcmb  = ft_getopt(cfg, 'channelcmb', {'all' 'all'});
cfg.refindx     = ft_getopt(cfg, 'refindx',    'all', 1);
cfg.trials      = ft_getopt(cfg, 'trials',     'all', 1);
cfg.complex     = ft_getopt(cfg, 'complex',    'abs');
cfg.jackknife   = ft_getopt(cfg, 'jackknife',  'no');
cfg.removemean  = ft_getopt(cfg, 'removemean', 'yes');
cfg.partchannel = ft_getopt(cfg, 'partchannel','');
cfg.parameter   = ft_getopt(cfg, 'parameter',  []);

hasjack = (isfield(data, 'method') && strcmp(data.method, 'jackknife')) || (isfield(data, 'dimord') && strcmp(data.dimord(1:6), 'rptjck'));
hasrpt  = (isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'rpt'))) || (isfield(data, 'avg') && isfield(data.avg, 'mom')) || (isfield(data, 'trial') && isfield(data.trial, 'mom')); % FIXME old-fashioned pcc data
dojack  = strcmp(cfg.jackknife, 'yes');
normrpt = 0; % default, has to be overruled e.g. in plv, because of single replicate normalisation
normpow = 1; % default, has to be overruled e.g. in csd

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  data = ft_selectdata(tmpcfg, data);
  [cfg, data] = rollback_provenance(cfg, data);
end

% select channels/channelcombination of interest and set the cfg-options accordingly
if isfield(data, 'label')
  selchan = cell(0, 1);
  if ~isempty(cfg.channelcmb) && ~isequal(cfg.channelcmb, {'all' 'all'}) && size(cfg.channelcmb,2)==2
    tmpcmb = ft_channelcombination(cfg.channelcmb, data.label);
    tmpchan = unique(tmpcmb(:));
    cfg.channelcmb = ft_channelcombination(cfg.channelcmb(:, 1:2), tmpchan, 1);
    selchan = [selchan; unique(cfg.channelcmb(:))];
  elseif ~isempty(cfg.channelcmb) && isequal(cfg.channelcmb, {'all' 'all'})
    cfg.channelcmb = ft_channelcombination(cfg.channelcmb, data.label, 1);
    selchan = [selchan; unique(cfg.channelcmb(:))];
  end
  
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  selchan = [selchan; cfg.channel];
  if ~isempty(cfg.partchannel)
    cfg.partchannel = ft_channelselection(cfg.partchannel, data.label);
    selchan = [selchan; cfg.partchannel];
  end
  tmpcfg = [];
  tmpcfg.channel = unique(selchan);
  data = ft_selectdata(tmpcfg, data);
  % restore the provenance information
  [cfg, data] = rollback_provenance(cfg, data);
elseif isfield(data, 'labelcmb')
  cfg.channel = ft_channelselection(cfg.channel, unique(data.labelcmb(:)));
  if ~isempty(cfg.partchannel)
    ft_error('partialization is only possible without linearly indexed bivariate data');
  end
  if ~isempty(cfg.channelcmb)
    % FIXME do something extra here
  end
  % FIXME call ft_selectdata
end

% FIXME check which methods require hasrpt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data bookkeeping - ensure that the input data is appropriate for the method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
needrpt = 1; % logical flag to specify whether (pseudo)-repetitions are required in the lower level connectivity function (can be singleton)
switch cfg.method
  case {'coh' 'csd'}
    if ~isempty(cfg.partchannel)
      if hasrpt && ~hasjack
        ft_warning('partialisation on single trial observations is not supported, removing trial dimension');
        try
          data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'}, 'cmbrepresentation', 'fullfast');
          inparam = 'crsspctrm';
          hasrpt = 0;
        catch
          ft_error('partial coherence/csd is only supported for input allowing for a all-to-all csd representation');
        end
      else
        %         try
        %           data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'}, 'cmbrepresentation', 'full');
        %           inparam = 'crsspctrm';
        %         catch
        %           ft_error('partial coherence/csd is only supported for input allowing for a all-to-all csd representation');
        %         end
        inparam = 'crsspctrm';
      end
    else
      data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq' 'source' 'source+mesh'});
      inparam = 'crsspctrm';
    end
    
    if strcmp(cfg.method, 'csd')
      normpow = 0;
      outparam = 'crsspctrm';
    elseif strcmp(cfg.method, 'coh')
      outparam = 'cohspctrm';
    end
    
    dtype = ft_datatype(data);
    switch dtype
      case 'source'
        if isempty(cfg.refindx), ft_error('indices of reference voxels need to be specified'); end
        % if numel(cfg.refindx)>1, ft_error('more than one reference voxel is not yet supported'); end
      otherwise
    end
    % FIXME think of accommodating partial coherence for source data with only a few references
  case {'wpli'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'wplispctrm';
    if hasjack, ft_error('to compute wpli, data should be in rpt format'); end
  case {'wpli_debiased'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'wpli_debiasedspctrm';
    if hasjack, ft_error('to compute wpli, data should be in rpt format'); end
  case {'ppc'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'ppcspctrm';
    if hasjack, ft_error('to compute ppc, data should be in rpt format'); end
  case {'wppc'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    outparam = 'wppcspctrm';
    if hasjack, ft_error('to compute wppc, data should be in rpt format'); end
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
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq' 'source' 'source+mesh'});
    dtype = ft_datatype(data);
    switch dtype
      case {'freq' 'freqmvar'}
        inparam = 'powcovspctrm';
      case {'source' 'source+mesh'}
        inparam = 'powcov';
        if isempty(cfg.refindx), ft_error('indices of reference voxels need to be specified'); end
        % if numel(cfg.refindx)>1, ft_error('more than one reference voxel is not yet supported'); end
      otherwise
    end
    outparam = [cfg.method, 'spctrm'];
  case {'granger' 'instantaneous_causality' 'total_interdependence' 'transfer' 'iis'}
    % create subcfg for the spectral factorization
    if ~isfield(cfg, 'granger')
      cfg.granger = [];
    end
    cfg.granger.conditional = ft_getopt(cfg.granger, 'conditional', 'no');
    cfg.granger.block       = ft_getopt(cfg.granger, 'block', []);
    cfg.granger.channelcmb  = ft_getopt(cfg.granger, 'channelcmb', cfg.channelcmb);
    cfg                     = removefields(cfg, 'channelcmb');
    data = ft_checkdata(data, 'datatype', {'mvar' 'freqmvar' 'freq'});
    inparam = {'transfer', 'noisecov', 'crsspctrm'};
    if strcmp(cfg.method, 'granger'),                 outparam = 'grangerspctrm'; end
    if strcmp(cfg.method, 'instantaneous_causality'), outparam = 'instantspctrm'; end
    if strcmp(cfg.method, 'total_interdependence'),   outparam = 'totispctrm';    end
    if strcmp(cfg.method, 'transfer'),                outparam = {'transfer' 'noisecov' 'crsspctrm'}; end
    if strcmp(cfg.method, 'iis'),                     outparam = 'iis'; end
    % check whether the frequency bins are more or less equidistant
    dfreq = diff(data.freq)./mean(diff(data.freq));
    assert(all(dfreq>0.999) && all(dfreq<1.001), ['non equidistant frequency bins are not supported for method ',cfg.method]);
    
  case {'ddtf'}
    data = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = {'transfer' 'crsspctrm'};
    outparam = [cfg.method, 'spctrm'];
  case {'dtf' 'pdc' 'gpdc'}
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
    inparam  = 'mom';
    outparam = 'powcorrspctrm';
  case {'mi' 'di' 'dfi'}
    % create the subcfg for the mutual information
    if ~isfield(cfg, cfg.method), cfg.(cfg.method) = []; end
    cfg.(cfg.method).method  = ft_getopt(cfg.(cfg.method), 'method',  'gcmi'); % default to the Gaussian Copula based method
    cfg.(cfg.method).numbin  = ft_getopt(cfg.(cfg.method), 'numbin',  10);
    cfg.(cfg.method).lags    = ft_getopt(cfg.(cfg.method), 'lags',    0);
    cfg.(cfg.method).montage = ft_getopt(cfg.(cfg.method), 'montage', []);
    cfg.(cfg.method).complex = ft_getopt(cfg.(cfg.method), 'complex', 'complex');
    cfg.(cfg.method).combinelags = ft_getopt(cfg.(cfg.method), 'combinelags', false);
    cfg.(cfg.method).feature     = ft_getopt(cfg.(cfg.method), 'feature',     []);
    cfg.(cfg.method).precondition = ft_getopt(cfg.(cfg.method), 'precondition', false);
    
    % what are the input requirements?
    data  = ft_checkdata(data, 'datatype', {'raw' 'timelock' 'freq' 'source'});
    dtype = ft_datatype(data);
    if strcmp(dtype, 'timelock')
      if ~isfield(data, 'trial')
        inparam = 'avg';
      else
        inparam = 'trial';
      end
      hasrpt = (isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'rpt')));
      
      cfg.refchannel = ft_getopt(cfg, 'refchannel', []);
      cfg.refindx    = ft_getopt(cfg, 'refindx',    []);
    elseif strcmp(dtype, 'raw')
      inparam = 'trial';
      hasrpt  = 1;
    
      cfg.refchannel = ft_getopt(cfg, 'refchannel', []);
      cfg.refindx    = ft_getopt(cfg, 'refindx',    []);
    elseif strcmp(dtype, 'freq')
      inparam = 'something';
    else
      inparam = 'something else';
    end
    outparam = cfg.method;
    needrpt  = 1;
  case 'laggedcoherence'
    data = ft_checkdata(data, 'datatype', {'freq'});
    if ~isfield(data, 'fourierspctrm')
      error('this connectivity method requires a ''fourierspctrm'' in the input data');
    end
    inparam  = 'lcrsspctrm';
    outparam = 'lcohspctrm';
    
    % create the subcfg for the lagged coherence
    if ~isfield(cfg, 'laggedcoherence'), cfg.laggedcoherence = []; end
    cfg.laggedcoherence.lags = ft_getopt(cfg.laggedcoherence, 'lags', []);
    cfg.laggedcoherence.timeresolved = false;
    
  otherwise
    ft_error('unknown method % s', cfg.method);
end

dtype = ft_datatype(data);

% FIXME throw an error if cfg.complex~='abs', and dojack==1
% FIXME throw an error if no replicates and cfg.method='plv'
% FIXME trial selection has to be implemented still

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data bookkeeping - check whether the required inparam is present in the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(~isfield(data, inparam)) || (isfield(data, 'crsspctrm') && (ischar(inparam) && strcmp(inparam, 'crsspctrm')))
  if iscell(inparam)
    % in the case of multiple inparams, use the first one to check the
    % input data (e.g. checking for 'transfer' for requested granger)
    inparam = inparam{1};
  end
  
  switch dtype
    case {'freq' 'freqmvar'}
      if strcmp(inparam, 'crsspctrm')
        if isfield(data, 'fourierspctrm')
          [data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm', 'crsspctrm', dtype, 'channelcmb', cfg.channelcmb, 'keeprpt', normrpt);
        elseif strcmp(inparam, 'crsspctrm') && isfield(data, 'powspctrm')
          % if input data is old-fashioned, i.e. contains powandcsd
          [data, powindx, hasrpt] = univariate2bivariate(data, 'powandcsd', 'crsspctrm', dtype, 'channelcmb', cfg.channelcmb, 'keeprpt', normrpt);
        elseif isfield(data, 'labelcmb')
          powindx = labelcmb2indx(data.labelcmb);
        else
          powindx = [];
        end
      elseif strcmp(inparam, 'lcrsspctrm')
        [data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm', 'lcrsspctrm', dtype, 'channelcmb', cfg.channelcmb, 'timeresolved', cfg.laggedcoherence.timeresolved, 'lags', cfg.laggedcoherence.lags);
      elseif strcmp(inparam, 'powcovspctrm')
        if isfield(data, 'powspctrm')
          [data, powindx] = univariate2bivariate(data, 'powspctrm', 'powcovspctrm', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'channelcmb', cfg.channelcmb, 'sqrtflag', strcmp(cfg.method, 'amplcorr'));
        elseif isfield(data, 'fourierspctrm')
          [data, powindx] = univariate2bivariate(data, 'fourierspctrm', 'powcovspctrm', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'channelcmb', cfg.channelcmb, 'sqrtflag', strcmp(cfg.method, 'amplcorr'));
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
        
        % convert the inparam back to cell-array in the case of granger
        if any(strcmp(cfg.method, {'granger' 'instantaneous_causality' 'total_interdependence' 'transfer' 'iis'}))
          inparam = {'transfer' 'noisecov' 'crsspctrm'};
          tmpcfg  = ft_checkconfig(cfg, 'createsubcfg', {'granger'});
          optarg  = ft_cfg2keyval(tmpcfg.granger);
        elseif strcmp(cfg.method, 'ddtf')
          inparam = {'transfer' 'crsspctrm'};
          tmpcfg  = ft_checkconfig(cfg, 'createsubcfg', {'ddtf'});
          optarg  = ft_cfg2keyval(tmpcfg.ddtf);
        else
          tmpcfg  = ft_checkconfig(cfg, 'createsubcfg', {cfg.method});
          optarg  = ft_cfg2keyval(tmpcfg.(cfg.method));
        end
        
        % compute the transfer matrix
        data   = ft_connectivity_csd2transfer(data, optarg{:});
      end
      
    case {'source' 'source+mesh'}
      if ischar(cfg.refindx) && strcmp(cfg.refindx, 'all')
        cfg.refindx = 1:size(data.pos,1);
      elseif ischar(cfg.refindx)
        ft_error('cfg.refindx should be a 1xN vector, or ''all''');
      end
      if strcmp(inparam, 'crsspctrm')
        [data, powindx, hasrpt] = univariate2bivariate(data, 'mom', 'crsspctrm', dtype, 'channelcmb', cfg.refindx, 'keeprpt', 0);
        % [data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm', 'crsspctrm', dtype, 0, cfg.refindx, [], 1);
      elseif strcmp(inparam, 'powcov')
        if isfield(data, 'pow')
          [data, powindx, hasrpt] = univariate2bivariate(data, 'pow', 'powcov', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'channelcmb', cfg.refindx, 'sqrtflag', strcmp(cfg.method, 'amplcorr'), 'keeprpt', 0);
        elseif isfield(data, 'mom')
          [data, powindx, hasrpt] = univariate2bivariate(data, 'mom', 'powcov', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'channelcmb', cfg.refindx, 'sqrtflag', strcmp(cfg.method, 'amplcorr'), 'keeprpt', 0);
        end
      end
      
    case 'comp'
      [data, powindx, hasrpt] = univariate2bivariate(data, 'trial', 'cov', dtype, 'demeanflag', strcmp(cfg.removemean, 'yes'), 'channelcmb', cfg.channelcmb, 'sqrtflag', false, 'keeprpt', 1);
      
  end % switch dtype
  
elseif (isfield(data, 'crsspctrm') && (ischar(inparam) && strcmp(inparam, 'crsspctrm')))
  % this means that there is a sparse crsspctrm in the data
  
else
  powindx = [];
end % ensure that the bivariate measure exists

% do some additional work if single trial normalisation is required
% for example when plv needs to be computed
if normrpt && hasrpt
  if strcmp(inparam, 'crsspctrm')
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
if ~isempty(cfg.partchannel) && (isfield(data, 'label') || isfield(data, 'labelcmb'))
  if isfield(data, 'label')
    label = data.label;
  elseif isfield(data, 'labelcmb')
    [indx, label] = labelcmb2indx(data.labelcmb);
  end
  allchannel = ft_channelselection(cfg.channel, label);
  pchanindx  = match_str(allchannel, cfg.partchannel);
  kchanindx  = setdiff(1:numel(allchannel), pchanindx);
  keepchn    = allchannel(kchanindx);
  
  cfg.pchanindx   = pchanindx;
  cfg.allchanindx = kchanindx;
  
  partstr = '';
  for k = 1:numel(cfg.partchannel)
    partstr = [partstr, '-', cfg.partchannel{k}];
  end
  for k = 1:numel(keepchn)
    keepchn{k} = [keepchn{k}, '\', partstr(2:end)];
  end
  if isfield(data, 'label')
    % update labels of the partialed channels
    data.label = keepchn;
  elseif isfield(data, 'labelcmb')
    for k = 1:numel(data.labelcmb)
      data.labelcmb{k} = [data.labelcmb{k}, '\', partstr(2:end)];
    end
  end
  
else
  cfg.pchanindx   = [];
  cfg.allchanindx = [];
end

% check if jackknife is required
if hasrpt && dojack && hasjack
  % do nothing
elseif hasrpt && dojack && ~ismember(cfg.method, {'wpli','wpli_debiased','ppc','wppc'})
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
elseif hasrpt && ~ismember(cfg.method, {'wpli','wpli_debiased','ppc','wppc','powcorr_ortho','mi','di','dfi'})% || needrpt)
  % create dof variable
  if isfield(data, 'dof')
    dof = data.dof;
  elseif isfield(data, 'cumtapcnt')
    dof = sum(data.cumtapcnt);
  end
  tmpcfg = [];
  tmpcfg.avgoverrpt = 'yes';
  tmpcfg.nanmean = 'yes';
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
    optarg = {'feedback', cfg.feedback, 'dojack', dojack, 'debias', strcmp(cfg.method, 'wpli_debiased')};
    [datout, varout, nrpt] = ft_connectivity_wpli(data.(inparam), optarg{:});
    
  case {'wppc' 'ppc'}
    % weighted pairwise phase consistency or pairwise phase consistency
    optarg = {'feedback', cfg.feedback, 'dojack', dojack, 'weighted', strcmp(cfg.method, 'wppc')};
    [datout, varout, nrpt] = ft_connectivity_ppc(data.(inparam), optarg{:});
    
  case 'plv'
    % phase locking value
    optarg = {'complex', cfg.complex, 'dimord', data.dimord, 'feedback', cfg.feedback, 'pownorm', normpow, 'hasjack', hasjack};
    if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'amplcorr'
    % amplitude correlation
    if isfield(data, 'dimord')
      dimord = data.dimord;
    else
      dimord = data.([inparam, 'dimord']);
    end
    optarg = {'feedback', cfg.feedback, 'dimord', dimord, 'complex', 'real', 'pownorm', 1, 'pchanindx', [], 'hasjack', hasjack};
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'powcorr'
    % power correlation
    if isfield(data, 'dimord')
      dimord = data.dimord;
    else
      dimord = data.([inparam, 'dimord']);
    end
    optarg = {'feedback', cfg.feedback, 'dimord', dimord, 'complex', 'real', 'pownorm', 1, 'pchanindx', [], 'hasjack', hasjack};
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'transfer'
    % the necessary stuff has already been computed
    datout   = data.transfer;
    noisecov = data.noisecov;
    crsspctrm = data.crsspctrm;
    if ~hasrpt
      datout   = shiftdim(datout,1);
      noisecov = shiftdim(noisecov,1);
      crsspctrm = shiftdim(crsspctrm,1);
    end
    
  case {'granger' 'instantaneous_causality' 'total_interdependence' 'iis'}
    % granger causality
    if ft_datatype(data, 'freq') || ft_datatype(data, 'freqmvar')
      if isfield(data, 'labelcmb') && isfield(cfg.granger, 'sfmethod') && strcmp(cfg.granger.sfmethod, 'bivariate_conditional')
        % create a powindx variable that ft_connectivity_granger can use to
        % do the conditioning
        [indx, label, blockindx, blocklabel] = labelcmb2indx(data.labelcmb);
        cmbindx12 = labelcmb2indx(cfg.granger.channelcmb(:,1:2), label);
        cmbindx23 = labelcmb2indx(cfg.granger.channelcmb(:,2:3), label);
        cmbindx   = [cmbindx12 cmbindx23(:,2); cmbindx12(:,[2 1]) cmbindx23(:,2)];
        
        powindx.cmbindx   = indx;
        powindx.blockindx = blockindx;
        powindx.outindx   = cmbindx;
        
        newlabelcmb = cell(size(cmbindx,1),2);
        for k = 1:size(newlabelcmb,1)
          newlabelcmb{k,1} = sprintf('%s|%s',label{cmbindx(k,2)},label{cmbindx(k,3)}); % deliberate swap of 2/1 as per the conventional definition in conditional granger computation
          newlabelcmb{k,2} = sprintf('%s|%s',label{cmbindx(k,1)},label{cmbindx(k,3)});
        end
        data.labelcmb = newlabelcmb;
        
      elseif isfield(data, 'labelcmb') && ~istrue(cfg.granger.conditional)
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
        
        if isfield(data, 'label')
          % this field should be removed
          data = rmfield(data, 'label');
        end
        
      elseif isfield(data, 'labelcmb') && istrue(cfg.granger.conditional)
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
        
        % make a temporary label list
        tmp2 = cell(numel(data.labelcmb),1);
        for m = 1:numel(data.labelcmb)
          tok = tokenize(data.labelcmb{m}, '[');
          tmp2{m} = tok{1};
        end
        label = cat(1,data.block.label); %unique(tmp2);
        
        [powindx.cmbindx, powindx.n] = blockindx2cmbindx(data.labelcmb, {label data.blockindx}, tmp);
        data.labelcmb                = newlabelcmb;
        
        if isfield(data, 'label')
          % this field should be removed
          data = rmfield(data, 'label');
        end
        
      elseif isfield(cfg.granger, 'block') && ~isempty(cfg.granger.block)
        % make a temporary label list
        if isfield(data, 'label')
          label = data.label;
        else
          tmp = cell(numel(data.labelcmb),1);
          for m = 1:numel(data.labelcmb)
            tok = tokenize(data.labelcmb{m}, '[');
            tmp{m} = tok{1};
          end
          label = unique(tmp);
        end
        
        % blockwise granger
        for k = 1:numel(cfg.granger.block)
          %newlabel{k, 1} = cat(2, cfg.granger.block(k).label{:});
          newlabel{k,1}  = cfg.granger.block(k).name;
          powindx{k,1}   = match_str(label, cfg.granger.block(k).label);
        end
        data.label = newlabel;
      else
        powindx = [];
      end
      if ~exist('powindx', 'var'), powindx = []; end
      if strcmp(cfg.method, 'granger'),                 methodstr = 'granger';      end
      if strcmp(cfg.method, 'instantaneous_causality'), methodstr = 'instantaneous'; end
      if strcmp(cfg.method, 'total_interdependence'),   methodstr = 'total';        end
      if strcmp(cfg.method, 'iis'),                     methodstr = 'iis';          end
      optarg = {'hasjack', hasjack, 'method', methodstr, 'powindx', powindx, 'dimord', data.dimord};
      [datout, varout, nrpt] = ft_connectivity_granger(data.transfer, data.noisecov, data.crsspctrm, optarg{:});
      if strcmp(cfg.method, 'iis')
        data.freq   = nan;
      end
    else
      ft_error('granger for time domain data is not yet implemented');
    end
    
  case {'dtf' 'ddtf'}
    % directed transfer function
    if isfield(data, 'labelcmb')
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
    optarg = {'feedback', cfg.feedback, 'powindx', powindx, 'hasjack', hasjack};
    hasrpt = ~isempty(strfind(data.dimord, 'rpt'));
    if hasrpt
      datin = data.transfer;
    else
      datin = reshape(data.transfer, [1 size(data.transfer)]);
      data.crsspctrm = reshape(data.crsspctrm, [1 size(data.crsspctrm)]);
    end
    if strcmp(cfg.method, 'ddtf'), optarg = cat(2, optarg, {'crsspctrm' data.crsspctrm}); end
    [datout, varout, nrpt] = ft_connectivity_dtf(datin, optarg{:});
    
  case {'pdc' 'gpdc'}
    % partial directed coherence
    if isfield(data, 'labelcmb')
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
    optarg = {'feedback', cfg.feedback, 'powindx', powindx, 'hasjack', hasjack};
    if strcmp(cfg.method, 'gpdc'), optarg = cat(2, optarg, {'noisecov' data.noisecov}); end
    hasrpt = ~isempty(strfind(data.dimord, 'rpt'));
    if hasrpt
      datin = data.(inparam);
    else
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
      posdim = find(strcmp(dimtok, '{pos}'));
      posdim = 4; % we concatenate across positions...
      rptdim = find(~cellfun('isempty',strfind(dimtok, 'rpt')));
      rptdim = rptdim-1; % the posdim has to be taken into account...
      dat    = cat(4, data.mom{data.inside});
      dat    = permute(dat,[posdim rptdim setdiff(1:ndims(dat),[posdim rptdim])]);
      
      datout = ft_connectivity_powcorr_ortho(dat, optarg{:});
      
      % HACK continued: format the output according to the inside and
      % refindx specifications
      if ischar(cfg.refindx) && strcmp(cfg.refindx, 'all')
        % create all-to-all output
        tmp = zeros(numel(data.inside));
        tmp(data.inside,data.inside) = datout;
        datout = tmp;
        clear tmp;
        
        outdimord = 'pos_pos_freq';
      else
        % create all-to-few output
        tmp = zeros(numel(data.inside), numel(cfg.refindx));
        tmp(data.inside, :) = datout;
        datout = tmp;
        clear tmp;
        
        outdimord = 'pos_pos_freq';
      end
    elseif strcmp(data.dimord, 'rpttap_chan_freq')
      % loop over all frequencies
      [nrpttap, nchan, nfreq] = size(data.fourierspctrm);
      datout = cell(1, nfreq);
      for i=1:length(data.freq)
        dat       = reshape(data.fourierspctrm(:,:,i), nrpttap, nchan).';
        datout{i} = ft_connectivity_powcorr_ortho(dat, optarg{:});
      end
      datout = cat(3, datout{:});
      % HACK otherwise I don't know how to inform the code further down about the dimord
      data.dimord = 'rpttap_chan_chan_freq';
    else
      ft_error('unsupported data representation');
    end
    varout = [];
    nrpt = numel(data.cumtapcnt);
    
  case {'mi' 'di' 'dfi'}
    % mutual information using the information breakdown toolbox, or gcmi
    % presence of the toolbox is checked in the low-level function.
    % directed information using the gcmi toolbox, requires a lag to be
    % specified
    if (strcmp(cfg.method, 'di') || strcmp(cfg.method, 'dfi')) && any(cfg.(cfg.method).lags<=0)
      error('directed information requires cfg.di.lags to be > 0');
    end
    
    if ~strcmp(dtype, 'raw') && (numel(cfg.(cfg.method).lags)>1 || cfg.(cfg.method).lags~=0)
      ft_error('computation of lagged mutual information is only possible with ''raw'' data in the input');
    end
    
    % if not specified prior to the call, make sure empty 'opts' field exists
    if ~isfield(cfg.(cfg.method), 'opts')
      cfg.(cfg.method).opts = [];
    end
 
    switch dtype
      case 'raw'
        % ensure the lags to be in samples, not in seconds.
        cfg.(cfg.method).lags = round(cfg.(cfg.method).lags.*data.fsample);
        
        % check which row(s) in the data are the reference
        if isempty(cfg.refchannel) && isempty(cfg.refindx)
          error('either ''cfg.refchannel'', or ''cfg.refindx'' should be specified');
        elseif ~isempty(cfg.refchannel)
          cfg.refindx = match_str(data.label, cfg.refchannel);
        elseif ischar(cfg.refindx) && strcmp(cfg.refindx, 'all')
          %error('this is yet not possible and should be fixed elegantly'); %FIXME now we should decide whether we allow for multivariate reference channels, or we treat the refindx as a per-element vector, i.e. allow for all-to-all
        end
        
        if strcmp(cfg.method, 'dfi') || strcmp(cfg.method, 'mi')
          cfg.(cfg.method).feature = ft_getopt(cfg.(cfg.method), 'feature', []);
          if strcmp(cfg.method, 'dfi') && isempty(cfg.dfi.feature)
            error('dfi requires a feature to be specified');
          end
          cfg.(cfg.method).featureindx = match_str(data.label, cfg.(cfg.method).feature);
          cfg.(cfg.method).featurelags = ft_getopt(cfg.(cfg.method), 'featurelags');
          if ~isempty(cfg.(cfg.method).featurelags), cfg.(cfg.method).featurelags = round(cfg.(cfg.method).featurelags.*data.fsample); end
        end
        
        dat = catnan(data.trial, max(abs(cfg.(cfg.method).lags)));
               
        % deal with cfg.mi.montage, which allows for multivariate stuff
        if ~isempty(cfg.(cfg.method).montage)
          [i1, i2]  = match_str(data.label, cfg.(cfg.method).montage.labelorg);
          i3        = setdiff((1:numel(data.label))',i1);
          tra       = cfg.(cfg.method).montage.tra(:,i2);
          tra(end+(1:numel(i3)),end+(1:numel(i3))) = eye(numel(i3));
          newlabel  = [cfg.(cfg.method).montage.labelnew;data.label(i3)];
          
          % update the refindx
          cfg.refindx = match_str(newlabel, cfg.refchannel);
          
          dat       = dat([i1; i3], :);
          refindx   = cfg.refindx; 
        else
          tra      = [];
          newlabel = [];
          refindx  = cfg.refindx;
        end
        
        if ischar(cfg.refindx) && strcmp(cfg.refindx, 'all')
          outdimord = 'chan_chan';
        elseif numel(cfg.refindx)==1 || numel(cfg.refchannel)==1,
          outdimord = 'chan';
        else
          %ft_error('at present cfg.refindx should be either ''all'', or scalar');
        end
        if numel(cfg.(cfg.method).lags)>1 && ~istrue(cfg.(cfg.method).combinelags)
          data.time = cfg.(cfg.method).lags./data.fsample;
          outdimord = [outdimord, '_time'];
        else
          data = rmfield(data, 'time');
        end
        
        if ~isempty(newlabel)
          data.label = newlabel;
        end
        
      case 'timelock'
        dat = data.(inparam);
        dat = reshape(permute(dat, [2 3 1]), [size(dat, 2) size(dat, 1)*size(dat, 3)]);
        data = rmfield(data, 'time');
        if ischar(cfg.refindx) && strcmp(cfg.refindx, 'all')
          outdimord = 'chan_chan';
        elseif numel(cfg.refindx)==1
          outdimord = 'chan';
        else
          ft_error('at present cfg.refindx should be either ''all'', or scalar');
        end
        
        %data.dimord = 'chan_chan';
      case 'freq'
        ft_error('not yet implemented');
        
      case 'source'
        % for the time being work with mom
        % dat = cat(2, data.mom{data.inside}).';
        dat = cat(1, data.mom{data.inside});
        % dat = abs(dat);
    end
    optarg = {'numbin', cfg.(cfg.method).numbin, 'lags',    cfg.(cfg.method).lags,    'refindx', refindx, ...
              'method', cfg.(cfg.method).method, 'complex', cfg.(cfg.method).complex, 'precondition', cfg.(cfg.method).precondition, ...
              'opts',   cfg.(cfg.method).opts};
    if ~isempty(tra),             optarg = cat(2, optarg, {'tra' tra});                                   end
    if strcmp(cfg.method, 'mi'),  optarg = cat(2, optarg, {'conditional', false});                        end
    if strcmp(cfg.method, 'mi'),  optarg = cat(2, optarg, {'featureindx', cfg.(cfg.method).featureindx}); end
    if strcmp(cfg.method, 'mi'),  optarg = cat(2, optarg, {'featurelags', cfg.(cfg.method).featurelags}); end
    if strcmp(cfg.method, 'mi'),  optarg = cat(2, optarg, {'combinelags', cfg.(cfg.method).combinelags}); end
    if strcmp(cfg.method, 'di'),  optarg = cat(2, optarg, {'conditional', true});                         end
    if strcmp(cfg.method, 'di'),  optarg = cat(2, optarg, {'combinelags', cfg.(cfg.method).combinelags}); end
    if strcmp(cfg.method, 'dfi'), optarg = cat(2, optarg, {'conditional', true});                         end
    if strcmp(cfg.method, 'dfi'), optarg = cat(2, optarg, {'featureindx', cfg.(cfg.method).featureindx}); end
    if strcmp(cfg.method, 'dfi'), optarg = cat(2, optarg, {'featurelags', cfg.(cfg.method).featurelags}); end
    if strcmp(cfg.method, 'dfi'), optarg = cat(2, optarg, {'combinelags', cfg.(cfg.method).combinelags}); end
    [datout] = ft_connectivity_mutualinformation(dat, optarg{:});
    varout   = [];
    nrpt     = [];    
  case 'corr'
    % pearson's correlation coefficient
    optarg = {'dimord', getdimord(data, inparam), 'feedback', cfg.feedback, 'hasjack', hasjack};
    if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'xcorr'
    % cross-correlation function
    ft_error('method %s is not yet implemented', cfg.method);
    
  case 'spearman'
    % spearman's rank correlation
    ft_error('method %s is not yet implemented', cfg.method);
    
  case 'laggedcoherence'
    % lagged coherence estimate
    optarg = {'complex', cfg.complex, 'dimord', data.dimord, 'feedback', cfg.feedback, 'pownorm', normpow, 'hasjack', hasjack};
    optarg = cat(2, optarg, {'powindx', powindx});
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    data = removefields(data, 'dof'); % the dof is not to be trusted
  
  otherwise
    ft_error('unknown method %s', cfg.method);
    
end % switch method

% remove the auto combinations if necessary -> FIXME this is granger specific and thus could move to ft_connectivity_granger
if (strcmp(cfg.method, 'granger') || strcmp(cfg.method, 'instantaneous_causality') || strcmp(cfg.method, 'total_interdependence')) && isfield(cfg, 'granger') && isfield(cfg.granger, 'sfmethod') && strcmp(cfg.granger.sfmethod, 'bivariate')
  % remove the auto-combinations based on the order in the data
  switch dtype
    case {'freq' 'freqmvar'}
      keepchn = 1:size(datout, 1);
      keepchn = mod(keepchn, 4)==2 | mod(keepchn, 4)==3;
      datout = datout(keepchn, :, :, :, :);
      if ~isempty(varout)
        varout = varout(keepchn, :, :, :, :);
      end
      data.labelcmb = data.labelcmb(keepchn, :);
    case 'source'
      % not yet implemented
  end
end

if exist('powindx', 'var') && ~isempty(powindx)
  % based on powindx
  switch dtype
    case {'freq' 'freqmvar'}
      if isfield(data, 'labelcmb') && ~isstruct(powindx)
        keepchn = powindx(:, 1) ~= powindx(:, 2);
        datout = datout(keepchn, :, :, :, :);
        if ~isempty(varout)
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
      if ~isempty(varout)
        varout = varout(keepchn, :, :, :, :);
      end
      inside = false(zeros(1, size(data.pos, 1)));
      inside(data.inside) = true;
      inside = inside(keepchn);
      %       data.inside = find(inside)';
      %       data.outside = find(inside==0)';
      data.pos = data.pos(keepchn, :);
  end % switch dtype
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the output structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch dtype
  case {'freq' 'freqmvar'}
    stat = keepfields(data, {'label', 'labelcmb', 'grad', 'elec', 'opto'});
    if isfield(data, 'labelcmb')
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
    if ~iscell(outparam)
      stat.(outparam) = datout;
      if ~isempty(varout)
        stat.([outparam, 'sem']) = (varout./nrpt).^0.5;
      end
    else
      stat.(outparam{1}) = datout;
      for k = 2:numel(outparam)
        if exist(outparam{k}, 'var')
          stat.(outparam{k}) = eval(outparam{k});
        end
      end
    end
    
  case 'timelock'
    stat = keepfields(data, {'label', 'labelcmb', 'grad', 'elec', 'opto'});
    % deal with the dimord
    if exist('outdimord', 'var')
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
    if ~isempty(varout)
      stat.([outparam, 'sem']) = (varout./nrpt).^0.5;
    end
    
  case {'source' 'source+mesh'}
    stat = keepfields(data, {'pos', 'dim', 'transform', 'inside', 'outside' 'tri'});
    stat.(outparam) = datout;
    if ~isempty(varout)
      stat.([outparam, 'sem']) = (varout/nrpt).^0.5;
    end
    
    % deal with the dimord
    if exist('outdimord', 'var')
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
    
  case 'raw'
    stat = [];
    stat.label = data.label;
    stat.(outparam) = datout;
    if ~isempty(varout)
      stat.([outparam, 'sem']) = (varout/nrpt).^0.5;
    end
    if exist('outdimord', 'var')
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
  if k==1
    begsmp = 1+nnans;
    endsmp = nsmp(1)+nnans;
  else
    begsmp = k*nnans + sum(nsmp(1:(k-1))) + 1;
    endsmp = k*nnans + sum(nsmp(1:k));
  end
  datamatrix(:,begsmp:endsmp) = datacells{k};
end
