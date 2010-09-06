function [stat] = ft_connectivityanalysis(cfg, data)

% FT_CONNECTIVITYANALYIS computes various measures of connectivity
% between MEG/EEG channels or between source-level signals.
%
% Use as
%   stat = ft_connectivityanalysis(cfg, data)
%   stat = ft_connectivityanalysis(cfg, timelock)
%   stat = ft_connectivityanalysis(cfg, freq)
%   stat = ft_connectivityanalysis(cfg, source)
% where the first input argument is a configuration structure (see
% below) and the second argument is the output of FT_PREPROCESSING,
% FT_TIMELOCKANLAYSIS or FT_FREQANALYSIS or FT_MVARANALYSIS or
% FT_SOURCEANALYSIS, depending on the connectivity metric that you
% want to compute.
%
% The configuration structure has to contain
%   cfg.method  = 'coh',       coherence, support for freq, freqmvar and
%                               source data. For partial coherence also
%                               specify cfg.partchannel
%                 'csd',       cross-spectral density matrix, can also
%                               calculate partial csds -
%                               if cfg.partchannel is specified
%                 'plv',       phase-locking value, support for freq and freqmvar data
%                 'corr',      correlation coefficient (Pearson)
%                 'xcorr',     cross correlation function
%                 'powcorr',   power correlation, support for freq and source data
%                 'amplcorr',  amplitude correlation, support for freq and source data
%                 'dtf',       directed transfer function, support for freq and freqmvar data
%                 'pdc',       partial directed coherence, support for freq and freqmvar data
%                 'granger',   granger causality, support for freq and freqmvar data
%                 'psi',       phaseslope index, support for freq and freqmvar data
%                 'pcd',       pairwise circular difference
%                 'di',        directionality index

% Copyright (C) 2009, Robert Oostenveld, Jan-Mathijs Schoffelen, Andre Bastos
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

%fieldtripdefs

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the defaults

%FIXME do method specific calls to checkconfig
if ~isfield(cfg, 'feedback'),   cfg.feedback   = 'none'; end
if ~isfield(cfg, 'channel'),    cfg.channel    = 'all'; end
if ~isfield(cfg, 'channelcmb'), cfg.channelcmb = {'all' 'all'};    end
if ~isfield(cfg, 'refindx'),    cfg.refindx    = [];    end
if ~isfield(cfg, 'trials'),     cfg.trials     = 'all'; end
if ~isfield(cfg, 'complex'),    cfg.complex    = 'abs'; end
if ~isfield(cfg, 'jackknife'),  cfg.jackknife  = 'no';  end
if ~isfield(cfg, 'removemean'), cfg.removemean = 'yes'; end
if ~isfield(cfg, 'partchannel'), cfg.partchannel = '';  end
if ~isfield(cfg, 'conditional'), cfg.conditional = [];  end
if ~isfield(cfg, 'blockindx'),   cfg.blockindx   = {};  end
if ~isfield(cfg, 'inputfile'),  cfg.inputfile = [];     end
if ~isfield(cfg, 'outputfile'), cfg.outputfile = [];    end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
  end
end

hasjack = (isfield(data, 'method') && strcmp(data.method, 'jackknife')) || (isfield(data, 'dimord') && strcmp(data.dimord(1:6), 'rptjck'));
hasrpt  = (isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'rpt'))) || (isfield(data, 'avg') && isfield(data.avg, 'mom')); %FIXME old-fashioned pcc data
dojack  = strcmp(cfg.jackknife, 'yes');
normrpt = 0; % default, has to be overruled e.g. in plv, because of single replicate normalisation
normpow = 1; % default, has to be overruled e.g. in csd,

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  data = selectdata(data, 'rpt', cfg.trials);
  if isfield(data, 'cfg'),
    cfg.trl    = findcfg(data.cfg, 'trl');
    cfg.trlold = findcfg(data.cfg, 'trlold');
  end
end

% FIXME check which methods require hasrpt

% ensure that the input data is appropriate for the method
switch cfg.method
  case {'coh' 'csd'}
    if ~isempty(cfg.partchannel)
      if hasrpt && ~hasjack,
        error('partialisation on single trial observations is not supported');
      end
      try
        data    = checkdata(data, 'datatype', {'freqmvar' 'freq'}, 'cmbrepresentation', 'full');
        inparam = 'crsspctrm';
      catch
        error('partial coherence/csd is only supported for input allowing for a all-to-all csd representation');
      end
    else
      data    = checkdata(data, 'datatype', {'freqmvar' 'freq' 'source'});
      inparam = 'crsspctrm';
    end
    
    if strcmp(cfg.method, 'csd'),
      normpow     = 0;
      % warning('cfg.complex for requested csd is set to %s, do you really want this?', cfg.complex);
    end
    
    dtype   = datatype(data);
    switch dtype
      case 'source'
        if isempty(cfg.refindx), error('indices of reference voxels need to be specified'); end
        % if numel(cfg.refindx)>1, error('more than one reference voxel is not yet supported'); end
      otherwise
    end
    % FIXME think of accommodating partial coherence for source data with only a few references
    
  case {'plv'}
    data    = checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
    normrpt = 1;
  case {'corr' 'xcorr'}
    data = checkdata(data, 'datatype', 'raw');
  case {'amplcorr' 'powcorr'}
    data    = checkdata(data, 'datatype', {'freqmvar' 'freq' 'source'});
    dtype   = datatype(data);
    switch dtype
      case {'freq' 'freqmvar'}
        inparam = 'powcovspctrm';
      case 'source'
        inparam = 'powcov';
        if isempty(cfg.refindx), error('indices of reference voxels need to be specified'); end
        % if numel(cfg.refindx)>1, error('more than one reference voxel is not yet supported'); end
      otherwise
    end
  case {'granger'}
    data    = checkdata(data, 'datatype', {'mvar' 'freqmvar' 'freq'});
    inparam = {'transfer', 'noisecov', 'crsspctrm'};
    % FIXME could also work with time domain data
  case {'instantaneous_causality'}
    data    = checkdata(data, 'datatype', {'mvar' 'freqmvar' 'freq'});
    inparam = {'transfer', 'noisecov', 'crsspctrm'};
  case {'total_interdependence'}
    data    = checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
  case {'dtf' 'pdc'}
    data    = checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'transfer';
  case {'psi'}
    if ~isfield(cfg, 'normalize'),  cfg.normalize  = 'no';  end
    data    = checkdata(data, 'datatype', {'freqmvar' 'freq'});
    inparam = 'crsspctrm';
  case {'di'}
    %wat eigenlijk?
  otherwise
    error('unknown method %s', cfg.method);
end
dtype = datatype(data);

% ensure that source data is in 'new' representation
if strcmp(dtype, 'source'),
  data = checkdata(data, 'sourcerepresentation', 'new');
end

% FIXME throw an error if cfg.complex~='abs', and dojack==1
% FIXME throw an error if no replicates and cfg.method='plv'
% FIXME trial selection has to be implemented still

if isfield(data, 'label'),
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  if ~isempty(cfg.partchannel)
    cfg.partchannel = ft_channelselection(cfg.partchannel, data.label);
  end
end

if isfield(data, 'label') && ~isempty(cfg.channelcmb),
  cfg.channelcmb = ft_channelcombination(cfg.channelcmb, cfg.channel, 1);
end

% check whether the required inparam is present in the data
if any(~isfield(data, inparam)) || (isfield(data, 'crsspctrm') && (ischar(inparam) && strcmp(inparam, 'crsspctrm'))),
  switch dtype
    case 'freq'
      if strcmp(inparam, 'crsspctrm')
        if isfield(data, 'fourierspctrm')
          [data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm', 'crsspctrm', dtype, 0, cfg.channelcmb);
        elseif strcmp(inparam, 'crsspctrm') && isfield(data, 'powspctrm')
          % if input data is old-fashioned, i.e. contains powandcsd
          [data, powindx, hasrpt] = univariate2bivariate(data, 'powandcsd', 'crsspctrm', dtype, 0, cfg.channelcmb);
        elseif isfield(data, 'labelcmb')
          powindx = labelcmb2indx(data.labelcmb);
        else 
          powindx = [];
        end
      elseif strcmp(inparam, 'powcovspctrm')
        if isfield(data, 'powspctrm'),
          [data, powindx] = univariate2bivariate(data, 'powspctrm', 'powcovspctrm', dtype, strcmp(cfg.removemean,'yes'), cfg.channelcmb, strcmp(cfg.method,'amplcorr'));
        elseif isfield(data, 'fourierspctrm'),
          [data, powindx] = univariate2bivariate(data, 'fourierspctrm', 'powcovspctrm', dtype, strcmp(cfg.removemean,'yes'), cfg.channelcmb, strcmp(cfg.method,'amplcorr'));
        end
      end
    case 'source'
      if strcmp(inparam, 'crsspctrm')
        [data, powindx, hasrpt] = univariate2bivariate(data, 'mom', 'crsspctrm', dtype, 0, cfg.refindx, [], 0);
        %[data, powindx, hasrpt] = univariate2bivariate(data, 'fourierspctrm', 'crsspctrm', dtype, 0, cfg.refindx, [], 1);
      elseif strcmp(inparam, 'powcov')
        data            = checkdata(data, 'haspow', 'yes');
        [data, powindx] = univariate2bivariate(data, 'pow', 'powcov', dtype, strcmp(cfg.removemean,'yes'), cfg.refindx, strcmp(cfg.method,'amplcorr'), 0);
      end
    otherwise
  end
  
else
  powindx = [];
end

% do some additional work if single trial normalisation is required
% for example when plv needs to be computed
if normrpt && hasrpt,
  if strcmp(inparam, 'crsspctrm'),
    tmp  = data.(inparam);
    nrpt = size(tmp,1);
    progress('init', cfg.feedback, 'normalising...');
    for k = 1:nrpt
      progress(k/nrpt, 'normalising amplitude of replicate %d from %d to 1\n', k, nrpt);
      tmp(k,:,:,:,:) = tmp(k,:,:,:,:)./abs(tmp(k,:,:,:,:));
    end
    progress('close');
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
elseif hasrpt && dojack,
  % compute leave-one-outs
  data    = selectdata(data, 'jackknife', 'yes');
  hasjack = 1;
elseif hasrpt
  data   = selectdata(data, 'avgoverrpt', 'yes');
  hasrpt = 0;
else
  % nothing required
end

% ensure that the first dimension is singleton if ~hasrpt
if hasrpt
  % nothing required
else
  if ischar(inparam)
    data.(inparam) = reshape(data.(inparam), [1 size(data.(inparam))]);
  else
    for k = 1:numel(inparam)
      data.(inparam{k}) = reshape(data.(inparam{k}), [1 size(data.(inparam{k}))]);
    end
  end
  data.dimord    = ['rpt_',data.dimord];
end
    
% compute the desired connectivity metric
switch cfg.method
  case 'coh'
    % coherence (unsquared), if cfg.complex = 'imag' imaginary part of
    % coherency
    
    tmpcfg             = [];
    tmpcfg.complex     = cfg.complex;
    tmpcfg.feedback    = cfg.feedback;
    tmpcfg.dimord      = data.dimord;
    tmpcfg.pownorm     = normpow;
    tmpcfg.pchanindx   = cfg.pchanindx;
    tmpcfg.allchanindx = cfg.allchanindx;
    tmpcfg.hasjack     = hasjack;
    if exist('powindx', 'var'), tmpcfg.powindx     = powindx; end
    optarg             = cfg2keyval(tmpcfg);
    
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg);
    outparam = 'cohspctrm';
    
  case 'csd'
    % cross-spectral density (only useful if partialisation is required)
    
    tmpcfg             = [];
    tmpcfg.complex     = cfg.complex;
    tmpcfg.feedback    = cfg.feedback;
    tmpcfg.dimord      = data.dimord;
    tmpcfg.pownorm     = normpow;
    tmpcfg.pchanindx   = cfg.pchanindx;
    tmpcfg.allchanindx = cfg.allchanindx;
    tmpcfg.hasjack     = hasjack;
    if exist('powindx', 'var'), tmpcfg.powindx     = powindx; end
    optarg             = cfg2keyval(tmpcfg);
    
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg);
    outparam = 'crsspctrm';
    
  case 'plv'
    % phase locking value
    
    tmpcfg           = [];
    tmpcfg.complex   = cfg.complex;
    tmpcfg.feedback  = cfg.feedback;
    tmpcfg.dimord    = data.dimord;
    tmpcfg.pownorm     = normpow;
    tmpcfg.pchanindx   = cfg.pchanindx;
    tmpcfg.allchanindx = cfg.allchanindx;
    tmpcfg.hasjack     = hasjack;
    if exist('powindx', 'var'), tmpcfg.powindx     = powindx; end
    optarg             = cfg2keyval(tmpcfg);
    
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg);
    outparam         = 'plvspctrm';
    
  case 'corr'
    % pearson's correlation coefficient
  case 'xcorr'
    % cross-correlation function
  case 'spearman'
    % spearman's rank correlation
  case 'amplcorr'
    % amplitude correlation
    
    tmpcfg             = [];
    tmpcfg.feedback    = cfg.feedback;
    tmpcfg.dimord      = data.dimord;
    tmpcfg.complex     = 'real';
    tmpcfg.pownorm     = 1;
    tmpcfg.pchanindx   = [];
    tmpcfg.hasjack     = hasjack;
    if exist('powindx', 'var'), tmpcfg.powindx = powindx; end
    optarg             = cfg2keyval(tmpcfg);
    
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg);
    outparam        = 'amplcorrspctrm';
    
  case 'powcorr'
    % power correlation
    
    tmpcfg             = [];
    tmpcfg.feedback    = cfg.feedback;
    tmpcfg.dimord      = data.dimord;
    tmpcfg.complex     = 'real';
    tmpcfg.hasjack     = hasjack;
    if exist('powindx', 'var'), tmpcfg.powindx = powindx; end
    optarg             = cfg2keyval(tmpcfg);
    
    [datout, varout, nrpt] = ft_connectivity_corr(data.(inparam), optarg);
    outparam        = 'powcorrspctrm';
    
  case 'granger'
    % granger causality
    
    if sum(datatype(data, {'freq' 'freqmvar'})),
    
      if isfield(data, 'labelcmb') && isempty(cfg.conditional),
        % multiple pairwise non-parametric transfer functions
        % linearly indexed
        powindx = labelcmb2indx(data.labelcmb);
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
        % do nothing
      end
      
      %fs = cfg.fsample; %FIXME do we really need this, or is this related to how
      %noisecov is defined and normalised?
      fs = 1;
      [datout, varout, nrpt] = ft_connectivity_granger(data.transfer, data.noisecov, data.crsspctrm, fs, hasjack, powindx);
      outparam = 'grangerspctrm';
    else
      error('granger for time domain data is not yet implemented');
    end
    
  case 'instantaneous_causality'
    % instantaneous ft_connectivity between the series, requires the same elements as granger
    if sum(datatype(data, {'freq' 'freqmvar'})),
      
      if isfield(data, 'labelcmb'),
        % multiple pairwise non-parametric transfer functions
        % linearly indexed
        powindx = labelcmb2indx(data.labelcmb);
      elseif isfield(cfg, 'block') && ~isempty(cfg.block)
        % blockwise granger
        powindx = cfg.block;
        for k = 1:2
          newlabel{k,1} = cat(2,powindx{k});
        end
        data.label = newlabel;
      else
        % do nothing
      end
      %fs = cfg.fsample; %FIXME do we really need this, or is this related to how
      %noisecov is defined and normalised?
      fs = 1;
      [datout, varout, nrpt] = ft_connectivity_instantaneous(data.transfer, data.noisecov, data.crsspctrm, fs, hasjack, powindx);
      outparam = 'instantspctrm';
    else
      error('instantaneous causality for time domain data is not yet implemented');
    end
    
  case 'total_interdependence'
    %total interdependence
    
    tmpcfg           = [];
    tmpcfg.complex   = 'abs';
    tmpcfg.feedback  = cfg.feedback;
    tmpcfg.dimord    = data.dimord;
    tmpcfg.pownorm   = normpow;
    tmpcfg.pchanindx = cfg.pchanindx;
    tmpcfg.allchanindx = cfg.allchanindx;
    tmpcfg.hasrpt      = hasrpt;
    tmpcfg.hasjack     = hasjack;
    if exist('powindx', 'var'), tmpcfg.powindx     = powindx; end
    optarg             = cfg2keyval(tmpcfg);
    
    [datout, varout, nrpt] = ft_connectivity_toti(tmpcfg, data.(inparam), hasrpt, hasjack);
    outparam         = 'totispctrm';
    
  case 'dtf'
    % directed transfer function
    
    if isfield(data, 'labelcmb'),
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
    
    tmpcfg          = [];
    tmpcfg.feedback = cfg.feedback;
    tmpcfg.powindx  = powindx;
    tmpcfg.hasjack  = hasjack;
    optarg          = cfg2keyval(tmpcfg);
    
    hasrpt = ~isempty(strfind(data.dimord, 'rpt'));
    if hasrpt,
      nrpt  = size(data.(inparam),1);
      datin = data.(inparam);
    else
      nrpt  = 1;
      datin = reshape(data.(inparam), [1 size(data.(inparam))]);
    end
    [datout, varout, nrpt] = ft_connectivity_dtf(datin, optarg);
    outparam = 'dtfspctrm';
    
  case 'pdc'
    % partial directed coherence
    
    if isfield(data, 'labelcmb'),
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
    
    tmpcfg          = [];
    tmpcfg.feedback = cfg.feedback;
    tmpcfg.powindx  = powindx;
    tmpcfg.hasjack  = hasjack;
    optarg          = cfg2keyval(tmpcfg);
    
    hasrpt = ~isempty(strfind(data.dimord, 'rpt'));
    if hasrpt,
      nrpt  = size(data.(inparam),1);
      datin = data.(inparam);
    else
      nrpt  = 1;
      datin = reshape(data.(inparam), [1 size(data.(inparam))]);
    end
    
    [datout, varout, nrpt] = ft_connectivity_pdc(datin, optarg);
    outparam = 'pdcspctrm';
    
  case 'pcd'
    % pairwise circular distance
  case 'psi'
    % phase slope index
    
    tmpcfg           = [];
    tmpcfg.feedback  = cfg.feedback;
    tmpcfg.dimord    = data.dimord;
    tmpcfg.nbin      = nearest(data.freq, data.freq(1)+cfg.bandwidth)-1;
    tmpcfg.normalize = cfg.normalize;
    tmpcfg.hasrpt      = hasrpt;
    tmpcfg.hasjack     = hasjack;
    if exist('powindx', 'var'), tmpcfg.powindx     = powindx; end
    optarg             = cfg2keyval(tmpcfg);
    
    [datout, varout, nrpt] = ft_connectivity_psi(data.(inparam), optarg);
    outparam         = 'psispctrm';
    
  case 'di'
    % directionality index
  otherwise
    error('unknown method %s', cfg.method);
end

%remove the auto combinations if necessary
if exist('powindx', 'var') && ~isempty(powindx),
  switch dtype
    case {'freq' 'freqmvar'}
      if isfield(data, 'labelcmb') && ~isstruct(powindx),
        keepchn = powindx(:,1) ~= powindx(:,2);
        datout  = datout(keepchn,:,:,:,:);
        if ~isempty(varout),
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
      inside = logical(zeros(1,size(data.pos,1)));
      inside(data.inside) = true;
      inside = inside(keepchn);
      data.inside  = find(inside)';
      data.outside = find(inside==0)';
      data.pos     = data.pos(keepchn,:);
  end
end

%create output structure
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
      stat.([outparam,'sem']) = (varout/nrpt).^0.5;
    end
  case 'source'
    stat         = [];
    stat.pos     = data.pos;
    stat.dim     = data.dim;
    stat.inside  = data.inside;
    stat.outside = data.outside;
    stat.(outparam) = datout;
    if ~isempty(varout),
      stat.([outparam,'sem']) = (varout/nrpt).^0.5;
    end
end

if isfield(data, 'freq'), stat.freq = data.freq; end
if isfield(data, 'frequency'), stat.frequency = data.frequency; end
if isfield(data, 'time'), stat.time = data.time; end
if isfield(data, 'grad'), stat.grad = data.grad; end
if isfield(data, 'elec'), stat.elec = data.elec; end
if exist('nrpt', 'var'),  stat.dof  = nrpt;      end
%FIXME this is not correct for TF-representations when trials have
%different lengths

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id   = '$Id$';

% remember the configuration details of the input data
try cfg.previous = data.cfg; end
% remember the exact configuration details in the output
stat.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', stat); % use the variable name "data" in the output file
end

%--------------------------------------------------------------
%function [c, v, n] = ft_connectivity_corr(cfg, input, hasrpt, hasjack)

%-------------------------------------------------------------
function [c, v, n] = ft_connectivity_toti(cfg, input, hasrpt, hasjack)

% FIXME move functionality into ft_connectivity_granger

cfg.hasrpt  = hasrpt;
cfg.hasjack = hasjack;
optarg = cfg2keyval(cfg);

[c, v, n] = ft_connectivity_corr(input, optarg);
c = -log(1-c.^2);
v = -log(1-v.^2); %FIXME this is probably not correct

%-------------------------------------------------------------
%function [c, v, n] = ft_connectivity_psi(cfg, input, hasrpt, hasjack)

%------------------------------------------------------------
%function [pdc, pdcvar, n] = ft_connectivity_pdc(cfg, input, hasjack)

%------------------------------------------------------------
%function [dtf, dtfvar, n] = ft_connectivity_dtf(cfg, input, hasjack)

%-------------------------------------------------------------------------
%function [granger, v, n] = ft_connectivity_granger(H, Z, S, fs, hasjack, powindx)

%----------------------------------------------------------------
function [instc, v, n] = ft_connectivity_instantaneous(H, Z, S, fs, hasjack,powindx)

% FIXME move functionality into ft_connectivity_granger

%Usage: causality = hz2causality(H,S,Z,fs);
%Inputs: transfer  = transfer function,
%        crsspctrm = 3-D spectral matrix;
%        noisecov  = noise covariance,
%        fs        = sampling rate
%Outputs: instantaneous causality spectrum between the channels.
%Total Interdependence = Granger (X->Y) + Granger (Y->X) + Instantaneous Causality
%               : auto-causality spectra are set to zero
% Reference: Brovelli, et. al., PNAS 101, 9849-9854 (2004), Rajagovindan
% and Ding, PLoS One Vol. 3, 11, 1-8 (2008)
%M. Dhamala, UF, August 2006.

%FIXME speed up code and check
siz = size(H);
if numel(siz)==4,
  siz(5) = 1;
end
n   = siz(1);
Nc  = siz(2);

outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));
if isempty(powindx)
  
  %clear S; for k = 1:size(H,3), h = squeeze(H(:,:,k)); S(:,:,k) = h*Z*h'/fs; end;
  for kk = 1:n
    for ii = 1:Nc
      for jj = 1:Nc
        if ii ~=jj,
          zc1     = reshape(Z(kk,jj,jj,:) - Z(kk,ii,jj,:).^2./Z(kk,ii,ii,:),[1 1 1 1 siz(5)]);
          zc1     = repmat(zc1,[1 1 1 siz(4) 1]);
          zc2     = reshape(Z(kk,ii,ii,:) - Z(kk,jj,ii,:).^2./Z(kk,jj,jj,:),[1 1 1 1 siz(5)]);
          zc2     = repmat(zc2,[1 1 1 siz(4) 1]);
          CTH1    = reshape(ctranspose(squeeze(H(kk,ii,jj,:,:))),1,1,1,siz(4));
          CTH2    = reshape(ctranspose(squeeze(H(kk,jj,ii,:,:))),1,1,1,siz(4));
          term1   = (S(kk,ii,ii,:,:) - H(kk,ii,jj,:,:).*zc1.*CTH1);
          term2   = (S(kk,jj,jj,:,:) - H(kk,jj,ii,:,:).*zc2.*CTH2);
          Sdet      = (S(kk,ii,ii,:,:).*S(kk,jj,jj,:,:)) - (S(kk,ii,jj,:,:).*S(kk,jj,ii,:,:));
          outsum(jj,ii,:) = outsum(jj,ii) + log((term1.*term2)./Sdet(kk,:,:,:));
          outssq(jj,ii,:) = outssq(jj,ii) + log((term1.*term2)./Sdet(kk,:,:,:)).^2;
        end
      end
      outsum(ii,ii,:,:) = 0;%self-granger set to zero
    end
  end
elseif ~iscell(powindx)
  % data are linearly indexed
  for k = 1:Nc
    for j = 1:n
      iauto1  = find(sum(powindx==powindx(k,1),2)==2);
      iauto2  = find(sum(powindx==powindx(k,2),2)==2);
      icross1 = k;
      icross2 = find(sum(powindx==powindx(ones(Nc,1)*k,[2 1]),2)==2);
      if iauto1 ~= iauto2
        zc1     = Z(j,iauto1) - Z(j,icross2).^2./Z(j,iauto2);
        zc1     = repmat(zc1,[1 1 siz(3)]);
        zc2     = Z(j,iauto2) - Z(j,icross1).^2./Z(j,iauto1);
        zc2     = repmat(zc2,[1 1 siz(3)]);
        CTH1    = reshape(ctranspose(squeeze(H(j,icross2,:))),1,1,siz(3));
        CTH2    = reshape(ctranspose(squeeze(H(j,icross1,:))),1,1,siz(3));
        term1   = (S(j,iauto2,:) - H(j,icross2,:).*zc1.*CTH1);
        term2   = (S(j,iauto1,:) - H(j,icross1,:).*zc2.*CTH2);
        Sdet      = (S(j,iauto2,:).*S(j,iauto1,:)) - (S(j,icross2,:).*S(j,icross1,:));
        outsum(k,:) = outsum(k) + log((term1.*term2)./Sdet(j,:,:));
        outssq(k,:) = outssq(k) + log((term1.*term2)./Sdet(j,:,:)).^2;
      end
    end
  end
  
  
end
instc = outsum./n;

if n>1,
  if hasjack
    bias = (n-1).^2;
  else
    bias = 1;
  end
  v = bias*(outssq - (outsum.^2)./n)./(n - 1);
else
  v = [];
end

%----------------------------------------
function [indx] = labelcmb2indx(labelcmb)

%identify the auto-combinations
ncmb = size(labelcmb,1);
indx = zeros(ncmb,2);

label = unique(labelcmb(:));
nchan = numel(label);
autoindx = zeros(nchan,1);
for k = 1:nchan
  sel1 = strcmp(label{k}, labelcmb(:,1));
  sel2 = strcmp(label{k}, labelcmb(:,2));
  autoindx = find(sel1 & sel2);
  
  indx(sel1,1) = autoindx;
  indx(sel2,2) = autoindx;
end


%------------------------------------------------------------------------------------------------------------------
function [data, powindx, hasrpt] = univariate2bivariate(data, inparam, outparam, dtype, demeanflag, cmb, sqrtflag, keeprpt)

if nargin<8, keeprpt    = 1;  end
if nargin<7, sqrtflag   = 0;  end
if nargin<6, cmb        = []; end
if nargin<5, demeanflag = 0;  end

switch dtype
  case 'freq'
    ncmb  = size(cmb,1);
    nchan = numel(data.label);
    getpowindx = 0;
    if ncmb==0,
      error('no channel combinations are specified');
    elseif ncmb==nchan^2 || ncmb==(nchan+1)*nchan*0.5,
      dofull = 1;
    else
      dofull = 0;
    end
    
    if strcmp(inparam, 'fourierspctrm') && strcmp(outparam, 'crsspctrm'),
      %fourier coefficients -> cross-spectral density
      if dofull
        data    = checkdata(data, 'cmbrepresentation', 'full');
      else
        data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
      end
    elseif strcmp(inparam, 'powandcsd') && strcmp(outparam, 'crsspctrm'),
      if ~isempty(cmb),
        data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
        
        %ensure getting powindx later on to prevent crash
        getpowindx = 1;
      else
        %data    = checkdata(data, 'cmbrepresentation', 'full');
        %this should not be possible
        error('cannot convert to a full csd representation');
        
        
      end
    elseif strcmp(inparam, 'fourierspctrm') && strcmp(outparam, 'powcovspctrm'),
      %fourier coefficients -> power covariance
      data = checkdata(data, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
      
      if sqrtflag, data.powspctrm = sqrt(data.powspctrm); end
      
      %get covariance by using checkdata
      if demeanflag,
        nrpt = size(data.powspctrm,1);
        mdat = nanmean(data.powspctrm,1);
        data.powspctrm = data.powspctrm - mdat(ones(1,nrpt),:,:,:,:,:);
      end
      data.fourierspctrm = data.powspctrm; %this is necessary for checkdata to work
      data.dimord        = ['rpttap',data.dimord(4:end)];
      data               = rmfield(data, 'powspctrm');
      data.cumtapcnt(:)  = 1;
      data.cumsumcnt(:)  = 1;
      if ncmb < (nchan-1)*nchan*0.5,
        data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
      else
        data    = checkdata(data, 'cmbrepresentation', 'full');
      end
      data.powcovspctrm = data.crsspctrm;
      data              = rmfield(data, 'crsspctrm');
    elseif strcmp(inparam, 'powspctrm') && strcmp(outparam, 'powcovspctrm'),
      %power-spectral density -> power covariance
      
      if sqrtflag, data.powspctrm = sqrt(data.powspctrm); end
      
      %get covariance by using checkdata
      if demeanflag,
        nrpt = size(data.powspctrm,1);
        mdat = nanmean(data.powspctrm,1);
        data.powspctrm = data.powspctrm - mdat(ones(1,nrpt),:,:,:,:,:);
      end
      data.fourierspctrm = data.powspctrm; %this is necessary for checkdata to work
      data.dimord        = ['rpttap',data.dimord(4:end)];
      data               = rmfield(data, 'powspctrm');
      data.cumtapcnt(:)  = 1;
      data.cumsumcnt(:)  = 1;
      if ncmb < (nchan-1)*nchan*0.5,
        data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
      else
        data    = checkdata(data, 'cmbrepresentation', 'full');
      end
      data.powcovspctrm = data.crsspctrm;
      data = rmfield(data, 'crsspctrm');
    else
      error('unknown conversion from univariate to bivariate representation');
    end
    
    if ~isempty(cmb) && (ncmb < (nchan-1)*nchan*0.5 || getpowindx==1),
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
  case 'source'
    ncmb = numel(cmb);
    
    if strcmp(inparam, 'pow') && strcmp(outparam, 'powcov'),
      [nrpt,nvox] = size(data.pow);
      if sqrtflag, data.pow = sqrt(data.pow); end
      if demeanflag,
        mdat = nanmean(data.pow,1);
        data.pow = data.pow - mdat(ones(1,nrpt),:); %FIXME only works for 1 frequency
      end
      
      data.powcov = [data.pow .* data.pow(:,ones(1,nvox)*cmb) data.pow.*data.pow];
      data        = rmfield(data, 'pow');
      powindx     = [nvox+(1:nvox) nvox+(1:nvox); cmb*ones(1,nvox) nvox+(1:nvox)]';
      
      data.pos    = [data.pos repmat(data.pos(cmb,:),[nvox 1]);data.pos data.pos];
      data.inside = [data.inside(:); data.inside(:)+nvox];
      data.outside = [data.outside(:); data.outside(:)+nvox];
      data.dim(2) = size(data.pos,1);
    elseif strcmp(inparam, 'mom') && strcmp(outparam, 'crsspctrm'),
      %get mom as rpttap_pos_freq matrix
      %FIXME this assumes only 1 freq bin
      mom = zeros(size(data.mom{data.inside(1)},1), size(data.pos,1));
      mom(:, data.inside) = cat(2, data.mom{data.inside});
      if keeprpt,
        [nrpt,nvox]    = size(mom);
        data.crsspctrm = [mom.*conj(mom(:,ones(1,nvox)*cmb)) abs(mom).^2];
        data           = rmfield(data, 'mom');
        data           = rmfield(data, 'momdimord');
        powindx     = [nvox+(1:nvox) nvox+(1:nvox); cmb*ones(1,nvox) nvox+(1:nvox)]';
        
        data.pos    = [data.pos repmat(data.pos(cmb,:),[nvox 1]);data.pos data.pos];
        data.inside = [data.inside(:); data.inside(:)+nvox];
        data.outside = [data.outside(:); data.outside(:)+nvox];
      elseif ncmb<size(mom,2)
        %do it computationally more efficient
        [nrpt,nvox]    = size(mom);
        
        data.crsspctrm = reshape((transpose(mom)*conj(mom(:,cmb)))./nrpt, [nvox*ncmb 1]);
        tmppow         = mean(abs(mom).^2)';
        data.crsspctrm = cat(1, data.crsspctrm, tmppow);
        tmpindx1       = transpose(ncmb*nvox + ones(ncmb+1,1)*(1:nvox));
        tmpindx2       = repmat(tmpindx1(cmb(:),end), [1 nvox])';
        tmpindx3       = repmat(cmb(:), [1 nvox])'; %expressed in original voxel indices
        powindx        = [tmpindx1(:) [tmpindx2(:);tmpindx1(:,end)]];
        
        data.pos       = [repmat(data.pos, [ncmb 1]) data.pos(tmpindx3(:),:); data.pos data.pos];
        data.inside    = data.inside(:)*ones(1,ncmb+1) + (ones(length(data.inside),1)*nvox)*(0:ncmb);
        data.inside    = data.inside(:);
        data.outside   = setdiff((1:nvox*(ncmb+1))', data.inside);
        if isfield(data, 'momdimord'),
          data.crsspctrmdimord = ['pos_',data.momdimord(14:end)];%FIXME this assumes dimord to be 'rpttap_...'
        end
        data           = rmfield(data, 'mom');
        data           = rmfield(data, 'momdimord');
      else
        [nrpt,nvox]    = size(mom);
        data.crsspctrm = (transpose(mom)*conj(mom))./nrpt;
        data           = rmfield(data, 'mom');
        data           = rmfield(data, 'momdimord');
        powindx        = [];
        data.crsspctrmdimord = 'pos_pos_freq'; %FIXME hard coded
      end
      data.dimord = data.crsspctrmdimord;
      clear mom;
    else
      error('unknown conversion from univariate to bivariate representation');
    end
  otherwise
end

hasrpt  = (isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'rpt')));
