function [stat] = ft_connectivityanalysis(cfg, data)

% FT_CONNECTIVITYANALYSIS computes various measures of connectivity
% between MEG/EEG channels or between source-level signals.
%
% Use as
%   stat = ft_connectivityanalysis(cfg, data)
%   stat = ft_connectivityanalysis(cfg, timelock)
%   stat = ft_connectivityanalysis(cfg, freq)
%   stat = ft_connectivityanalysis(cfg, source)
% where the first input argument is a configuration structure (see
% below) and the second argument is the output of FT_PREPROCESSING,
% FT_TIMELOCKANLAYSIS, FT_FREQANALYSIS, FT_MVARANALYSIS,
% FT_SOURCEANALYSIS. The different connectivity
% metrics are supported only for specific datatypes (see below). 
%
% The configuration structure has to contain
%   cfg.method  =  string, can be
%     'coh',       coherence, support for freq, freqmvar and source data.
%                  For partial coherence also specify cfg.partchannel
%     'csd',       cross-spectral density matrix, can also calculate partial
%                  csds - if cfg.partchannel is specified, support for freq
%                  and freqmvar data
%     'plv',       phase-locking value, support for freq and freqmvar data
%     'powcorr',   power correlation, support for freq and source data
%     'amplcorr',  amplitude correlation, support for freq and source data
%     'granger',   granger causality, support for freq and freqmvar data
%     'dtf',       directed transfer function, support for freq and
%                  freqmvar data 
%     'pdc',       partial directed coherence, support for freq and 
%                  freqmvar data
%     'psi',       phaseslope index, support for freq and freqmvar data
%     'wpli',      weighted phase lag index (signed one,
%                  still have to take absolute value to get indication of
%                  strength of interaction. Note: measure has positive
%                  bias. Use wpli_debiased to avoid this.
%     'wpli_debiased'  debiased weighted phase lag index
%                  (estimates squared wpli)
%     'ppc'        pairwise phase consistency 
%     'wppc'       weighted pairwise phase consistency
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
% See also FT_PREPROCESSING, FT_TIMELOCKANALYSIS, FT_FREQANALYSIS,
% FT_MVARANALYSIS, FT_SOURCEANALYSIS, FT_NETWORKANALYSIS

% Undocumented options:
%   cfg.refindx     
%   cfg.conditional  
%   cfg.blockindx    
%   cfg.jackknife    

% Methods to be implemented
%
%                 'xcorr',     cross correlation function
%                 'di',        directionality index
%                 'spearman'   spearman's rank correlation
%                 'corr'       pearson correlation
%
% Copyright (C) 2009, Jan-Mathijs Schoffelen, Andre Bastos, Martin Vinck, Robert Oostenveld
% Copyright (C) 2010-2011, Jan-Mathijs Schoffelen, Martin Vinck
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

% FIXME it should be checked carefully whether the following works
% data = ft_checkdata(data, 'datatype', {'raw', 'timelock', 'freq', 'source'});

% set the defaults
cfg.feedback    = ft_getopt(cfg, 'feedback',    'none');
cfg.channel     = ft_getopt(cfg, 'channel',     'all');
cfg.channelcmb  = ft_getopt(cfg, 'channelcmb',  {'all' 'all'});
cfg.refindx     = ft_getopt(cfg, 'refindx',     'all');
cfg.trials      = ft_getopt(cfg, 'trials',      'all');
cfg.complex     = ft_getopt(cfg, 'complex',     'abs');
cfg.jackknife   = ft_getopt(cfg, 'jackknife',   'no');
cfg.removemean  = ft_getopt(cfg, 'removemean',  'yes');
cfg.partchannel = ft_getopt(cfg, 'partchannel', '');
cfg.inputfile   = ft_getopt(cfg, 'inputfile',   []);
cfg.outputfile  = ft_getopt(cfg, 'outputfile',  []);
cfg.parameter   = ft_getopt(cfg, 'parameter',   []);

hasjack = (isfield(data, 'method') && strcmp(data.method, 'jackknife')) || (isfield(data, 'dimord') && strcmp(data.dimord(1:6), 'rptjck'));
hasrpt  = (isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'rpt'))) || ...
  (isfield(data, 'avg') && isfield(data.avg, 'mom')); %FIXME old-fashioned pcc data
dojack  = strcmp(cfg.jackknife, 'yes');
normrpt = 0; % default, has to be overruled e.g. in plv, because of single replicate normalisation
normpow = 1; % default, has to be overruled e.g. in csd,

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  data = ft_selectdata(data, 'rpt', cfg.trials);
end

% select channels/channelcombination of interest and set the cfg-options accordingly
if isfield(data, 'label'),
  selchan = cell(0,1);
  if ~isempty(cfg.channelcmb) && ~isequal(cfg.channelcmb, {'all' 'all'}),
    tmpcmb         = ft_channelcombination(cfg.channelcmb, data.label);
    tmpchan        = unique(tmpcmb(:));
    cfg.channelcmb = ft_channelcombination(cfg.channelcmb, tmpchan, 1);
    selchan        = [selchan;unique(cfg.channelcmb(:))];
  end
  
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  selchan     = [selchan;cfg.channel];
  if ~isempty(cfg.partchannel)
    cfg.partchannel = ft_channelselection(cfg.partchannel, data.label);
    selchan         = [selchan; cfg.partchannel];
  end
  data = ft_selectdata(data, 'channel', unique(selchan));
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data bookkeeping:
% ensure that the input data is appropriate for the method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
needrpt = 1; % logical flag to specify whether (pseudo)-repetitions are required in the lower level
% connectivity function (can be singleton)
switch cfg.method
  case {'coh' 'csd'}
    if ~isempty(cfg.partchannel)
      if hasrpt && ~hasjack,
        warning('partialisation on single trial observations is not supported, removing trial dimension');
        try
          data    = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'}, 'cmbrepresentation', 'fullfast');
          inparam = 'crsspctrm';
          hasrpt  = 0;
        catch
          error('partial coherence/csd is only supported for input allowing for a all-to-all csd representation');
        end
      else
        try
          data    = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'}, 'cmbrepresentation', 'full');
          inparam = 'crsspctrm';
        catch
          error('partial coherence/csd is only supported for input allowing for a all-to-all csd representation');
        end
      end
    else
      data    = ft_checkdata(data, 'datatype', {'freqmvar' 'freq' 'source'});
      inparam = 'crsspctrm';
    end
    
    if strcmp(cfg.method, 'csd'),
      normpow     = 0;
      outparam    = 'crsspctrm';
    elseif strcmp(cfg.method, 'coh'),
      outparam    = 'cohspctrm';
    end
    
    dtype   = ft_datatype(data);
    switch dtype
      case 'source'
        if isempty(cfg.refindx), error('indices of reference voxels need to be specified'); end
        % if numel(cfg.refindx)>1, error('more than one reference voxel is not yet supported'); end
      otherwise
    end
    % FIXME think of accommodating partial coherence for source data with only a few references
  case {'wpli'}
    data     = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam  = 'crsspctrm';
    outparam = 'wplispctrm';
    debiaswpli = 0;
    if hasjack, error('to compute wpli, data should be in rpt format'); end
  case {'wpli_debiased'}
    data     = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam  = 'crsspctrm';
    outparam = 'wpli_debiasedspctrm';
    debiaswpli = 1;
    if hasjack, error('to compute wpli, data should be in rpt format'); end
  case {'ppc'}
    data     = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam  = 'crsspctrm';
    outparam = 'ppcspctrm';
    weightppc = 0;
    if hasjack, error('to compute ppc, data should be in rpt format'); end
  case {'wppc'}
    data     = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam  = 'crsspctrm';
    outparam = 'wppcspctrm';
    weightppc = 1;
    if hasjack, error('to compute wppc, data should be in rpt format'); end
  case {'plv'}
    data     = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam  = 'crsspctrm';
    outparam = 'plvspctrm';
    normrpt  = 1;
  case {'corr' 'xcorr'}
    data = ft_checkdata(data, 'datatype', 'raw');
  case {'amplcorr' 'powcorr'}
    data    = ft_checkdata(data, 'datatype', {'freqmvar' 'freq' 'source'});
    dtype   = ft_datatype(data);
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
  case {'granger'}
    if ~isfield(cfg, 'granger')
      cfg.granger = [];
    end
    cfg.granger.conditional = ft_getopt(cfg.granger, 'conditional', []);
    cfg.granger.blockindx   = ft_getopt(cfg.granger, 'blockindx',   {});
    if isfield(cfg, 'channelcmb'),
      cfg.granger.channelcmb = cfg.channelcmb;
      cfg = rmfield(cfg, 'channelcmb');
    end
    data     = ft_checkdata(data, 'datatype', {'mvar' 'freqmvar' 'freq'});
    inparam  = {'transfer', 'noisecov', 'crsspctrm'};
    outparam = 'grangerspctrm';
    % FIXME could also work with time domain data
  case {'instantaneous_causality'}
    data     = ft_checkdata(data, 'datatype', {'mvar' 'freqmvar' 'freq'});
    inparam  = {'transfer', 'noisecov', 'crsspctrm'};
    outparam = 'instantspctrm';
  case {'total_interdependence'}
    data     = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam  = {'transfer', 'noisecov', 'crsspctrm'};
    outparam = 'totispctrm';
  case {'dtf' 'pdc'}
    data     = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam  = 'transfer';
    outparam = [cfg.method,'spctrm'];
  case {'psi'}
    if ~isfield(cfg, 'normalize'),  cfg.normalize  = 'no';  end
    data     = ft_checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam  = 'crsspctrm';
    outparam = 'psispctrm';
  case {'hipp'}
    data     = ft_checkdata(data, 'datatype', 'source');
    %inparam  = 'avg.mom';
    inparam  = 'mom';
    outparam = 'powcorrspctrm';
  case {'di'}
    %wat eigenlijk?
  otherwise
    error('unknown method %s', cfg.method);
end
dtype = ft_datatype(data);

% ensure that source data is in 'new' representation
if strcmp(dtype, 'source'),
  data = ft_checkdata(data, 'sourcerepresentation', 'new');
end

% FIXME throw an error if cfg.complex~='abs', and dojack==1
% FIXME throw an error if no replicates and cfg.method='plv'
% FIXME trial selection has to be implemented still

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data bookkeeping:
% check whether the required inparam is present in the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(~isfield(data, inparam)) || (isfield(data, 'crsspctrm') && (ischar(inparam) && strcmp(inparam, 'crsspctrm'))),
  if iscell(inparam)
    % in the case of multiple inparams, use the first one to check the
    % input data (e.g. checking for 'transfer' for requested granger)
    inparam    = inparam{1};
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
          [data, powindx] = univariate2bivariate(data, 'powspctrm', 'powcovspctrm', dtype, 'demeanflag', strcmp(cfg.removemean,'yes'), 'cmb', cfg.channelcmb, 'sqrtflag', strcmp(cfg.method,'amplcorr'));
        elseif isfield(data, 'fourierspctrm'),
          [data, powindx] = univariate2bivariate(data, 'fourierspctrm', 'powcovspctrm', dtype, 'demeanflag', strcmp(cfg.removemean,'yes'), 'cmb', cfg.channelcmb, 'sqrtflag', strcmp(cfg.method,'amplcorr'));
        end
      elseif strcmp(inparam, 'transfer')
        if isfield(data, 'fourierspctrm')
          %FIXME this is fast but throws away the trial dimension, consider
          %a way to keep trial information if needed, but use the fast way
          %if possible
          data   = ft_checkdata(data, 'cmbrepresentation', 'fullfast');
          hasrpt = 0;
        elseif isfield(data, 'powspctrm')
          data = ft_checkdata(data, 'cmbrepresentation', 'full');
        end
        
        tmpcfg = ft_checkconfig(cfg, 'createsubcfg',  {'granger'});
        
        % check whether multiple pairwise decomposition is required (this
        % can most conveniently be handled at this level
        %tmpcfg.npsf = rmfield(tmpcfg.npsf, 'channelcmb');
        try,tmpcfg.npsf = rmfield(tmpcfg.granger, 'block');     end
        try,tmpcfg.npsf = rmfield(tmpcfg.granger, 'blockindx'); end
        %         end
        optarg = ft_cfg2keyval(tmpcfg.granger);
        data   = ft_connectivity_csd2transfer(data, optarg{:});
        
        % convert the inparam back to cell array in the case of granger
        if strcmp(cfg.method, 'granger') || strcmp(cfg.method, 'instantaneous_causality') || strcmp(cfg.method, 'total_interdependence')
          inparam = {'transfer' 'noisecov' 'crsspctrm'};
        end
      end
    case 'source'
      if strcmp(inparam, 'crsspctrm')
        [data, powindx, hasrpt] = univariate2bivariate(data, 'mom', 'crsspctrm', dtype, 'cmb', cfg.refindx, 'keeprpt', 0);
        %[data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm', 'crsspctrm', dtype, 0, cfg.refindx, [], 1);
      elseif strcmp(inparam, 'powcov')
        data            = ft_checkdata(data, 'sourcerepresentation', 'new', 'haspow', 'yes');
        [data, powindx, hasrpt] = univariate2bivariate(data, 'pow', 'powcov', dtype, 'demeanflag', strcmp(cfg.removemean,'yes'), 'cmb', cfg.refindx, 'sqrtflag', strcmp(cfg.method,'amplcorr'), 'keeprpt', 0);
      end
    otherwise
  end
  
elseif (isfield(data, 'crsspctrm') && (ischar(inparam) && strcmp(inparam, 'crsspctrm')))
  % this means that there is a sparse crsspctrm in the data
else
  powindx = [];
end

% do some additional work if single trial normalisation is required
% for example when plv needs to be computed
if normrpt && hasrpt,
  if strcmp(inparam, 'crsspctrm'),
    tmp  = data.(inparam);
    nrpt = size(tmp,1);
    ft_progress('init', cfg.feedback, 'normalising...');
    for k = 1:nrpt
      ft_progress(k/nrpt, 'normalising amplitude of replicate %d from %d to 1\n', k, nrpt);
      tmp(k,:,:,:,:) = tmp(k,:,:,:,:)./abs(tmp(k,:,:,:,:));
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
  pchanindx  = match_str(allchannel,cfg.partchannel);
  kchanindx  = setdiff(1:numel(allchannel), pchanindx);
  keepchn    = allchannel(kchanindx);
  
  cfg.pchanindx   = pchanindx;
  cfg.allchanindx = kchanindx;
  partstr = '';
  for k = 1:numel(cfg.partchannel)
    partstr = [partstr,'-',cfg.partchannel{k}];
  end
  for k = 1:numel(keepchn)
    keepchn{k} = [keepchn{k},'\',partstr(2:end)];
  end
  data.label      = keepchn; % update labels to remove the partialed channels
  % FIXME consider keeping track of which channels have been partialised
else
  cfg.pchanindx   = [];
  cfg.allchanindx = [];
end

% check if jackknife is required
if hasrpt && dojack && hasjack,
  % do nothing
elseif hasrpt && dojack && ~(exist('debiaswpli', 'var') || exist('weightppc', 'var')),
  % compute leave-one-outs
  data    = ft_selectdata(data, 'jackknife', 'yes');
  hasjack = 1;
elseif hasrpt && ~(exist('debiaswpli', 'var') || exist('weightppc', 'var') || strcmp(cfg.method, 'hipp'))
  % create dof variable
  if isfield(data, 'dof')
    dof = data.dof;
  elseif isfield(data, 'cumtapcnt')
    dof = sum(data.cumtapcnt);
  end
  data   = ft_selectdata(data, 'avgoverrpt', 'yes');
  hasrpt = 0;
else
  % nothing required
end

% ensure that the first dimension is singleton if ~hasrpt
if hasrpt
  % nothing required
elseif needrpt
  if ischar(inparam)
    data.(inparam) = reshape(data.(inparam), [1 size(data.(inparam))]);
  else
    for k = 1:numel(inparam)
      data.(inparam{k}) = reshape(data.(inparam{k}), [1 size(data.(inparam{k}))]);
    end
  end
  
  if isfield(data, 'dimord')
    data.dimord    = ['rpt_',data.dimord];
  else
    data.([inparam,'dimord']) = ['rpt_',data.([inparam,'dimord'])];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the desired connectivity metric by
% calling the appropriate ft_connectivity_XXX
% function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch cfg.method
  case 'coh'
    % coherence (unsquared), if cfg.complex = 'imag' imaginary part of coherency
    optarg = {'complex',  cfg.complex, 'dimord',  data.dimord, 'feedback', cfg.feedback, ...
      'pownorm',  normpow,     'hasjack', hasjack};
    if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'csd'
    % cross-spectral density (e.g. useful if partialisation is required)
    optarg = {'complex',  cfg.complex, 'dimord',  data.dimord, 'feedback', cfg.feedback, ...
      'pownorm',  normpow,     'hasjack', hasjack};
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
    optarg = {'complex',  cfg.complex, 'dimord',  data.dimord, 'feedback', cfg.feedback, ...
      'pownorm',  normpow,     'hasjack', hasjack};
    if ~isempty(cfg.pchanindx), optarg = cat(2, optarg, {'pchanindx', cfg.pchanindx, 'allchanindx', cfg.allchanindx}); end
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'corr'
    % pearson's correlation coefficient
  case 'xcorr'
    % cross-correlation function
  case 'spearman'
    % spearman's rank correlation
  case 'amplcorr'
    % amplitude correlation
    
    if isfield(data, 'dimord'),
      dimord = data.dimord;
    else
      dimord = data.([inparam,'dimord']);
    end
    optarg = {'feedback', cfg.feedback, 'dimord',    dimord, 'complex', 'real', ...
      'pownorm',  1,            'pchanindx', [],     'hasjack', hasjack};
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'powcorr'
    % power correlation
    
    if isfield(data, 'dimord'),
      dimord = data.dimord;
    else
      dimord = data.([inparam,'dimord']);
    end
    optarg = {'feedback', cfg.feedback, 'dimord',    dimord, 'complex', 'real', ...
      'pownorm',  1,            'pchanindx', [],     'hasjack', hasjack};
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg{:});
    
  case 'granger'
    % granger causality
    
    if sum(ft_datatype(data, {'freq' 'freqmvar'})),
      
      if isfield(data, 'labelcmb') && isempty(cfg.granger.conditional),
        % multiple pairwise non-parametric transfer functions
        % linearly indexed
        
        % The following is very slow, one may make assumptions regarding
        % the order of the channels -> csd2transfer gives combinations in
        % quadruplets, where the first and fourth are auto-combinations,
        % and the second and third are cross-combinations
        % powindx = labelcmb2indx(data.labelcmb);
        
        % The following is not needed anymore, because ft_connectivity_granger
        % relies on some hard-coded assumptions for the channel-pair ordering.
        % Otherwise it becomes just too slow.
        %         powindx = zeros(size(data.labelcmb));
        %         for k = 1:size(powindx,1)/4
        %           ix = ((k-1)*4+1):k*4;
        %           powindx(ix,:) = [1 1;4 1;1 4;4 4] + (k-1)*4;
        %         end
        
        powindx = [];
        
      elseif isfield(data, 'labelcmb')
        % conditional (blockwise) needs linearly represented cross-spectra
        for k = 1:size(cfg.conditional,1)
          tmp{k,1} = cfg.conditional(k,:);
          tmp{k,2} = cfg.conditional(k,[1 3]);
        end
        [cmbindx, n] = blockindx2cmbindx(data.labelcmb, cfg.blockindx, tmp);
        powindx.cmbindx = cmbindx;
        powindx.n       = n;
      elseif isfield(cfg, 'block') && ~isempty(cfg.block)
        % blockwise granger
        powindx = cfg.block;
        for k = 1:2
          newlabel{k,1} = cat(2,powindx{k});
        end
        data.label = newlabel;
      else
        powindx = [];
      end
      
      %fs = cfg.fsample; %FIXME do we really need this, or is this related to how
      %noisecov is defined and normalised?
      if ~exist('powindx', 'var'), powindx = []; end
      optarg = {'hasjack', hasjack, 'method', 'granger', 'powindx', powindx, 'dimord', data.dimord};
      [datout, varout, nrpt] = ft_connectivity_granger(data.transfer, data.noisecov, data.crsspctrm, optarg{:});
    else
      error('granger for time domain data is not yet implemented');
    end
    
  case 'instantaneous_causality'
    % instantaneous ft_connectivity between the series, requires the same elements as granger
    
    if sum(ft_datatype(data, {'freq' 'freqmvar'})),
      
      if isfield(data, 'labelcmb') && isempty(cfg.conditional),
        % linearly indexed channel pairs
      elseif isfield(data, 'labelcmb')
        % conditional (blockwise) needs linearly represented cross-spectra
        for k = 1:size(cfg.conditional,1)
          tmp{k,1} = cfg.conditional(k,:);
          tmp{k,2} = cfg.conditional(k,[1 3]);
        end
        [cmbindx, n] = blockindx2cmbindx(data.labelcmb, cfg.blockindx, tmp);
        powindx.cmbindx = cmbindx;
        powindx.n       = n;
      elseif isfield(cfg, 'block') && ~isempty(cfg.block)
        % blockwise granger
        powindx  = cfg.block;
        newlabel = cell(2,1);
        for k = 1:2
          newlabel{k,1} = cat(2,powindx{k});
        end
        data.label = newlabel;
      else
        powindx = [];
      end
      
      %fs = cfg.fsample; %FIXME do we really need this, or is this related to how
      %noisecov is defined and normalised?
      if ~exist('powindx', 'var'), powindx = []; end
      optarg = {'hasjack', hasjack, 'method', 'instantaneous', 'powindx', powindx, 'dimord', data.dimord};
      [datout, varout, nrpt] = ft_connectivity_granger(data.transfer, data.noisecov, data.crsspctrm, optarg{:});
    else
      error('granger for time domain data is not yet implemented');
    end
    
  case 'total_interdependence'
    %total interdependence
    if sum(ft_datatype(data, {'freq' 'freqmvar'})),
      
      if isfield(data, 'labelcmb') && isempty(cfg.conditional),
        % multiple pairwise non-parametric transfer functions
        % linearly indexed
      elseif isfield(data, 'labelcmb')
        % conditional (blockwise) needs linearly represented cross-spectra
        for k = 1:size(cfg.conditional,1)
          tmp{k,1} = cfg.conditional(k,:);
          tmp{k,2} = cfg.conditional(k,[1 3]);
        end
        [cmbindx, n] = blockindx2cmbindx(data.labelcmb, cfg.blockindx, tmp);
        powindx.cmbindx = cmbindx;
        powindx.n       = n;
      elseif isfield(cfg, 'block') && ~isempty(cfg.block)
        % blockwise granger
        powindx = cfg.block;
        for k = 1:2
          newlabel{k,1} = cat(2,powindx{k});
        end
        data.label = newlabel;
      else
        powindx = [];
      end
      
      %fs = cfg.fsample; %FIXME do we really need this, or is this related to how
      %noisecov is defined and normalised?
      if ~exist('powindx', 'var'), powindx = []; end
      optarg = {'hasjack', hasjack, 'method', 'total', 'powindx', powindx, 'dimord', data.dimord};
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
      nrpt  = size(data.(inparam),1);
      datin = data.(inparam);
    else
      nrpt  = 1;
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
      nrpt  = size(data.(inparam),1);
      datin = data.(inparam);
    else
      nrpt  = 1;
      datin = reshape(data.(inparam), [1 size(data.(inparam))]);
    end
    
    [datout, varout, nrpt] = ft_connectivity_pdc(datin, optarg{:});
    
  case 'psi'
    % phase slope index
    
    nbin   = nearest(data.freq, data.freq(1)+cfg.bandwidth)-1;
    optarg = {'feedback',  cfg.feedback,  'dimord', data.dimord, 'nbin',    nbin, ...
      'normalize', cfg.normalize, 'hasrpt', hasrpt,      'hasjack', hasjack};
    if exist('powindx', 'var'), optarg = cat(2, optarg, {'powindx', powindx}); end
    [datout, varout, nrpt] = ft_connectivity_psi(data.(inparam), optarg{:});
    
  case 'hipp'
    % Joerg Hipp's power correlation method
    
    optarg   = {'refindx', cfg.refindx, 'tapvec', data.cumtapcnt};
    [datout] = ft_connectivity_hipp(cat(2,data.mom{data.inside}).', optarg{:});
    varout   = [];
    nrpt     = numel(data.cumtapcnt);
    
  case 'di'
    % directionality index
  otherwise
    error('unknown method %s', cfg.method);
end

%remove the auto combinations if necessary -> FIXME this is granger specific and
%thus could move to ft_connectivity_granger
if (strcmp(cfg.method, 'granger') || strcmp(cfg.method, 'instantaneous_causality') || strcmp(cfg.method, 'total_interdependence')) && isfield(cfg, 'granger') && isfield(cfg.granger, 'sfmethod') && strcmp(cfg.granger.sfmethod, 'bivariate'),
  % remove the auto-combinations based on the order in the data
  switch dtype
    case {'freq' 'freqmvar'}
      keepchn = 1:size(datout,1);
      keepchn = mod(keepchn,4)==2 | mod(keepchn,4)==3;
      datout  = datout(keepchn,:,:,:,:);
      if ~isempty(varout),
        varout = varout(keepchn,:,:,:,:);
      end
      data.labelcmb = data.labelcmb(keepchn,:);
    case 'source'
      % not yet implemented
  end
end

if exist('powindx', 'var') && ~isempty(powindx),
  % based on powindx
  switch dtype
    case {'freq' 'freqmvar'}
      if isfield(data, 'labelcmb') && ~isstruct(powindx),
        keepchn = powindx(:,1) ~= powindx(:,2);
        datout  = datout(keepchn,:,:,:,:);
        if ~isempty(varout),
          if all(size(varout)==size(nrpt))
            nrpt = nrpt(keepchn,:,:,:,:);
          end
          varout = varout(keepchn,:,:,:,:);
        end
        data.labelcmb = data.labelcmb(keepchn,:);
      end
    case 'source'
      nvox    = size(unique(data.pos(:,1:3),'rows'),1);
      ncmb    = size(data.pos,1)/nvox-1;
      remove  = (powindx(:,1) == powindx(:,2)) & ((1:size(powindx,1))' > nvox*ncmb);
      keepchn = ~remove;
      
      datout = datout(keepchn,:,:,:,:);
      if ~isempty(varout),
        varout = varout(keepchn,:,:,:,:);
      end
      inside = false(zeros(1,size(data.pos,1)));
      inside(data.inside) = true;
      inside = inside(keepchn);
      data.inside  = find(inside)';
      data.outside = find(inside==0)';
      data.pos     = data.pos(keepchn,:);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the output structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch dtype
  case {'freq' 'freqmvar'},
    stat        = [];
    if isfield(data, 'label'),
      stat.label  = data.label;
    end
    if isfield(data, 'labelcmb'),
      stat.labelcmb = data.labelcmb;
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
      stat.([outparam,'sem']) = (varout./nrpt).^0.5;
    end
  case 'source'
    stat         = [];
    stat.pos     = data.pos;
    if isfield(stat, 'dim'),
      stat.dim     = data.dim;
    end
    stat.inside  = data.inside;
    stat.outside = data.outside;
    stat.(outparam) = datout;
    if ~isempty(varout),
      stat.([outparam,'sem']) = (varout/nrpt).^0.5;
    end
end

if isfield(data, 'freq'), stat.freq = data.freq; end
if isfield(data, 'time'), stat.time = data.time; end
if isfield(data, 'grad'), stat.grad = data.grad; end
if isfield(data, 'elec'), stat.elec = data.elec; end
if exist('dof',  'var'),  stat.dof  = dof;       end
% FIXME this needs to be implemented still

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data
ft_postamble history stat
ft_postamble savevar stat
