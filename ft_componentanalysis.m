function [comp] = ft_componentanalysis(cfg, data)

% FT_COMPONENTANALYSIS performs independent component analysis or other
% spatio-temporal decompositions of EEG or MEG data. This function computes
% the topography and timecourses of the components. The output of this
% function can be further analyzed with FT_TIMELOCKANALYSIS or
% FT_FREQANALYSIS.
%
% Use as
%   [comp] = ft_componentanalysis(cfg, data)
% where cfg is a configuration structure and the input data is obtained from
% FT_PREPROCESSING or from FT_TIMELOCKANALYSIS.
%
% The configuration should contain
%   cfg.method       = 'runica', 'fastica', 'binica', 'pca', 'svd', 'jader', 'varimax', 'dss', 'cca', 'sobi', 'white' or 'csp' (default = 'runica')
%   cfg.channel      = cell-array with channel selection (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.trials       = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.numcomponent = 'all' or number (default = 'all')
%   cfg.demean       = 'no' or 'yes', whether to demean the input data (default = 'yes')
%   cfg.updatesens   = 'no' or 'yes' (default = 'yes')
%   cfg.feedback     = 'no', 'text', 'textbar', 'gui' (default = 'text')
%
% The runica method supports the following method-specific options. The values that
% these options can take can be found with HELP RUNICA.
%   cfg.runica.extended
%   cfg.runica.pca
%   cfg.runica.sphering
%   cfg.runica.weights
%   cfg.runica.lrate
%   cfg.runica.block
%   cfg.runica.anneal
%   cfg.runica.annealdeg
%   cfg.runica.stop
%   cfg.runica.maxsteps
%   cfg.runica.bias
%   cfg.runica.momentum
%   cfg.runica.specgram
%   cfg.runica.posact
%   cfg.runica.verbose
%   cfg.runica.logfile
%   cfg.runica.interput
%
% The fastica method supports the following method-specific options. The values that
% these options can take can be found with HELP FASTICA.
%   cfg.fastica.approach
%   cfg.fastica.numOfIC
%   cfg.fastica.g
%   cfg.fastica.finetune
%   cfg.fastica.a1
%   cfg.fastica.a2
%   cfg.fastica.mu
%   cfg.fastica.stabilization
%   cfg.fastica.epsilon
%   cfg.fastica.maxNumIterations
%   cfg.fastica.maxFinetune
%   cfg.fastica.sampleSize
%   cfg.fastica.initGuess
%   cfg.fastica.verbose
%   cfg.fastica.displayMode
%   cfg.fastica.displayInterval
%   cfg.fastica.firstEig
%   cfg.fastica.lastEig
%   cfg.fastica.interactivePCA
%   cfg.fastica.pcaE
%   cfg.fastica.pcaD
%   cfg.fastica.whiteSig
%   cfg.fastica.whiteMat
%   cfg.fastica.dewhiteMat
%   cfg.fastica.only
%
% The binica method supports the following method-specific options. The values that
% these options can take can be found with HELP BINICA.
%   cfg.binica.extended
%   cfg.binica.pca
%   cfg.binica.sphering
%   cfg.binica.lrate
%   cfg.binica.blocksize
%   cfg.binica.maxsteps
%   cfg.binica.stop
%   cfg.binica.weightsin
%   cfg.binica.verbose
%   cfg.binica.filenum
%   cfg.binica.posact
%   cfg.binica.annealstep
%   cfg.binica.annealdeg
%   cfg.binica.bias
%   cfg.binica.momentum
%
% The dss method requires the following method-specific option and supports
% a whole lot of other options. The values that these options can take can be
% found with HELP DSS_CREATE_STATE.
%   cfg.dss.denf.function
%   cfg.dss.denf.params
%
% The sobi method supports the following method-specific options. The values that
% these options can take can be found with HELP SOBI.
%   cfg.sobi.n_sources
%   cfg.sobi.p_correlations
%
% The csp method implements the common-spatial patterns method. For CSP, the
% following specific options can be defined:
%   cfg.csp.classlabels = vector that assigns a trial to class 1 or 2.
%   cfg.csp.numfilters  = the number of spatial filters to use (default: 6).
%
% The icasso method implements icasso. It runs fastica a specified number of
% times, and provides information about the stability of the components found
% The following specific options can be defined, see ICASSOEST:
%   cfg.icasso.mode
%   cfg.icasso.Niter
%
% Instead of specifying a component analysis method, you can also specify
% a previously computed unmixing matrix, which will be used to estimate the
% component timecourses in this data. This requires
%   cfg.unmixing     = NxN unmixing matrix
%   cfg.topolabel    = Nx1 cell-array with the channel labels
%
% You may specify a particular seed for random numbers called by
% rand/randn/randi, or the random state used by a previous call to this
% function to replicate results. For example:
%   cfg.randomseed   = integer seed value of user's choice
%   cfg.randomseed   = comp.cfg.callinfo.randomseed (from previous call)
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_TOPOPLOTIC, FT_REJECTCOMPONENT, FASTICA, RUNICA, BINICA, SVD,
% JADER, VARIMAX, DSS, CCA, SOBI, ICASSO

% Copyright (C) 2003-2012, Robert Oostenveld
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

% undocumented cfg options:
%   cfg.cellmode = string, 'no' or 'yes', allows to run in cell-mode, i.e.
%     no concatenation across trials is needed. This is based on experimental
%     code and only supported for 'dss', 'fastica' and 'bsscca' as methods.

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
ft_preamble randomseed

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
istimelock = ft_datatype(data, 'timelock');
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'detrend'});
cfg = ft_checkconfig(cfg, 'renamed',    {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'predetermined mixing matrix', 'predetermined unmixing matrix'});
cfg = ft_checkconfig(cfg, 'deprecated', {'topo'});

% set the defaults
cfg.method          = ft_getopt(cfg, 'method',       'runica');
cfg.demean          = ft_getopt(cfg, 'demean',       'yes');
cfg.trials          = ft_getopt(cfg, 'trials',       'all', 1);
cfg.channel         = ft_getopt(cfg, 'channel',      'all');
cfg.numcomponent    = ft_getopt(cfg, 'numcomponent', 'all');
cfg.normalisesphere = ft_getopt(cfg, 'normalisesphere', 'yes');
cfg.cellmode        = ft_getopt(cfg, 'cellmode',     'no');
cfg.doscale         = ft_getopt(cfg, 'doscale',      'yes');
cfg.updatesens      = ft_getopt(cfg, 'updatesens',   'yes');
cfg.feedback        = ft_getopt(cfg, 'feedback',     'text');

% select channels, has to be done prior to handling of previous (un)mixing matrix
cfg.channel = ft_channelselection(cfg.channel, data.label);

if istrue(cfg.cellmode)
  ft_hastoolbox('cellfunction', 1);
end

if isfield(cfg, 'topo') && isfield(cfg, 'topolabel')
  ft_warning(['Specifying cfg.topo (= mixing matrix) to determine component '...
    'timecourses in specified data is deprecated; please specify an '...
    'unmixing matrix instead with cfg.unmixing. '...
    'Using cfg.unmixing=pinv(cfg.topo) for now to reproduce old behaviour.']);
  
  cfg.unmixing = pinv(cfg.topo);
  cfg = rmfield(cfg, 'topo');
end

if isfield(cfg, 'unmixing') && isfield(cfg, 'topolabel')
  % use the previously determined unmixing matrix on this dataset
  
  % test whether all required channels are present in the data
  [datsel, toposel] = match_str(cfg.channel, cfg.topolabel);
  if length(toposel)~=length(cfg.topolabel)
    ft_error('not all channels that are required for the unmixing are present in the data');
  end
  
  % ensure that all data channels not used in the unmixing should be removed from the channel selection
  tmpchan = match_str(cfg.channel, cfg.topolabel);
  cfg.channel = cfg.channel(tmpchan);
  
  % update some settings where there is no further choice to be made by the user
  cfg.numcomponent = 'all';
  cfg.method       = 'predetermined unmixing matrix';
end

% add the options for the specified methods to the configuration, only if needed
switch cfg.method
  case 'icasso'
    cfg.icasso        = ft_getopt(cfg,        'icasso', []);
    cfg.icasso.mode   = ft_getopt(cfg.icasso, 'mode',   'both');
    cfg.icasso.Niter  = ft_getopt(cfg.icasso, 'Niter',  15);
    cfg.icasso.method = ft_getopt(cfg.icasso, 'method', 'fastica');
    
    cfg.fastica       = ft_getopt(cfg, 'fastica', []);
  case 'fastica'
    % additional options, see FASTICA for details
    cfg.fastica = ft_getopt(cfg, 'fastica', []);
  case 'runica'
    % additional options, see RUNICA for details
    cfg.runica       = ft_getopt(cfg,        'runica',  []);
    cfg.runica.lrate = ft_getopt(cfg.runica, 'lrate',   0.001);
  case 'binica'
    % additional options, see BINICA for details
    cfg.binica       = ft_getopt(cfg,        'binica',  []);
    cfg.binica.lrate = ft_getopt(cfg.binica, 'lrate',   0.001);
  case {'dss' 'dss2'} % JM at present has his own dss, that can deal with cell-array input, specify as dds2
    % additional options, see DSS for details
    cfg.dss               = ft_getopt(cfg,          'dss',      []);
    cfg.dss.denf          = ft_getopt(cfg.dss,      'denf',     []);
    cfg.dss.denf.function = ft_getopt(cfg.dss.denf, 'function', 'denoise_fica_tanh');
    cfg.dss.denf.params   = ft_getopt(cfg.dss.denf, 'params',   []);
    cfg.dss.preprocf      = ft_getopt(cfg.dss,      'preprocf', []);
    cfg.dss.preprocf.function = ft_getopt(cfg.dss.preprocf, 'function', 'pre_sphere');
    cfg.dss.preprocf.params   = ft_getopt(cfg.dss.preprocf, 'params', []);
  case 'csp'
    % additional options, see CSP for details
    cfg.csp = ft_getopt(cfg, 'csp', []);
    cfg.csp.numfilters = ft_getopt(cfg.csp, 'numfilters', 6);
    cfg.csp.classlabels = ft_getopt(cfg.csp, 'classlabels');
  case 'bsscca'
    % additional options, see BSSCCA for details
    cfg.bsscca           = ft_getopt(cfg,        'bsscca', []);
    cfg.bsscca.refdelay  = ft_getopt(cfg.bsscca, 'refdelay', 1);
    cfg.bsscca.chandelay = ft_getopt(cfg.bsscca, 'chandelay', 0);
    if strcmp(cfg.cellmode, 'no')
      ft_error('cfg.mehod = ''bsscca'' requires cfg.cellmode = ''yes''');
    end
  otherwise
    % do nothing
end

% select trials of interest
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

Ntrials  = length(data.trial);
Nchans   = length(data.label);
if Nchans==0
  ft_error('no channels were selected');
end

% default is to compute just as many components as there are channels in the data
if strcmp(cfg.numcomponent, 'all')
  defaultNumCompsUsed = true(1);
  cfg.numcomponent = length(data.label);
else
  defaultNumCompsUsed = false(1);
end

% determine the size of each trial, they can be variable length
Nsamples = zeros(1,Ntrials);
for trial=1:Ntrials
  Nsamples(trial) = size(data.trial{trial},2);
end

if strcmp(cfg.demean, 'yes')
  % optionally perform baseline correction on each trial
  ft_info('baseline correcting data \n');
  for trial=1:Ntrials
    data.trial{trial} = ft_preproc_baselinecorrect(data.trial{trial});
  end
end

if strcmp(cfg.doscale, 'yes')
  % determine the scaling of the data, scale it to approximately unity
  % this will improve the performance of some methods, esp. fastica
  tmp                 = data.trial{1};
  tmp(~isfinite(tmp)) = 0; % ensure that the scaling is a finite value
  scale = norm((tmp*tmp')./size(tmp,2)); clear tmp;
  scale = sqrt(scale);
  if scale ~= 0
    ft_info('scaling data with 1 over %f\n', scale);
    for trial=1:Ntrials
      data.trial{trial} = data.trial{trial} ./ scale;
    end
  else
    ft_info('no scaling applied, since factor is 0\n');
  end
else
  ft_info('no scaling applied to the data\n');
end

if strcmp(cfg.method, 'sobi')
  
  % concatenate all the data into a 3D matrix respectively 2D (sobi)
  ft_info('concatenating data');
  Nsamples = Nsamples(1);
  dat = zeros(Ntrials, Nchans, Nsamples);
  % all trials should have an equal number of samples
  % and it is assumed that the time axes of all trials are aligned
  for trial=1:Ntrials
    ft_info('.');
    dat(trial,:,:) = data.trial{trial};
  end
  ft_info('\n');
  ft_info('concatenated data matrix size %dx%dx%d\n', size(dat,1), size(dat,2), size(dat,3));
  if Ntrials == 1
    dummy = 0;
    [dat, dummy] = shiftdim(dat);
  else
    dat = shiftdim(dat,1);
  end
  
elseif strcmp(cfg.method, 'csp')
  
  % concatenate the trials into two data matrices, one for each class
  sel1 = find(cfg.csp.classlabels==1);
  sel2 = find(cfg.csp.classlabels==2);
  if min(length(sel1), length(sel2)) == 0
    ft_error('CSP requires class labels!');
  end
  if length(sel1)+length(sel2)~=length(cfg.csp.classlabels)
    ft_warning('not all trials belong to class 1 or 2');
  end
  dat1 = cat(2, data.trial{sel1});
  dat2 = cat(2, data.trial{sel2});
  ft_info('concatenated data matrix size for class 1 is %dx%d\n', size(dat1,1), size(dat1,2));
  ft_info('concatenated data matrix size for class 2 is %dx%d\n', size(dat2,1), size(dat2,2));
  
elseif ~strcmp(cfg.method, 'predetermined unmixing matrix') && strcmp(cfg.cellmode, 'no')
  % concatenate all the data into a 2D matrix unless we already have an
  % unmixing matrix or unless the user request it otherwise
  ft_info('concatenating data');
  
  dat = zeros(Nchans, sum(Nsamples));
  for trial=1:Ntrials
    ft_info('.');
    begsample = sum(Nsamples(1:(trial-1))) + 1;
    endsample = sum(Nsamples(1:trial));
    dat(:,begsample:endsample) = data.trial{trial};
  end
  ft_info('\n');
  ft_info('concatenated data matrix size %dx%d\n', size(dat,1), size(dat,2));
  
  hasdatanans = any(~isfinite(dat(:)));
  if hasdatanans
    ft_info('data contains nans, only using the non-nan samples\n');
    finitevals = sum(~isfinite(dat))==0;
    dat        = dat(:,finitevals);
  end
else
  ft_info('not concatenating data\n');
  dat = data.trial;
  % FIXME cellmode processing is not nan-transparent yet
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the component analysis
ft_info('starting decomposition using %s\n', cfg.method);
switch cfg.method
  
  case 'icasso'
    % check whether the required low-level toolboxes are installed
    ft_hastoolbox('icasso',  1);
    
    if strcmp(cfg.icasso.method, 'fastica')
      ft_hastoolbox('fastica', 1);
      cfg.fastica.numOfIC = cfg.numcomponent;
      
      optarg = ft_cfg2keyval(cfg.(cfg.icasso.method));
      sR     = icassoEst(cfg.icasso.mode, dat, cfg.icasso.Niter, optarg{:});
    elseif strcmp(cfg.icasso.method, 'dss')
      % recurse into ft_componentanalysis and do some post processing
      tmpcfg        = rmfield(cfg, 'icasso');
      tmpcfg.method = cfg.icasso.method;
      tmpdata       = data;
      
      % initialize the variables to hold the output
      sR.W     = cell(cfg.icasso.Niter, 1);
      sR.A     = cell(cfg.icasso.Niter, 1);
      sR.index = zeros(0,2);
      for k = 1:cfg.icasso.Niter
        tmp = ft_componentanalysis(tmpcfg, tmpdata);
        sR.W{k}  = tmp.unmixing;
        sR.A{k}  = tmp.topo;
        sR.index = cat(1, sR.index, [k*ones(size(tmp.topo,2),1) (1:size(tmp.topo,2))']);
        sR.whiteningMatrix   = tmp.cfg.dss.V;
        sR.dewhiteningMatrix = tmp.cfg.dss.dV;
      end
      sR.signal = dat;
      sR.mode   = cfg.icasso.mode;
      sR.rdim   = size(tmp.topo,2);
    else
      ft_error('only ''fastica'' or ''dss'' is supported as method for icasso');
    end
    
    % do the rest of the icasso related processing
    sR = icassoCluster(sR, 'strategy', 'AL', 'simfcn', 'abscorr', 's2d', 'sim2dis', 'L',cfg.numcomponent);
    sR = icassoProjection(sR, 'cca', 's2d', 'sqrtsim2dis', 'epochs', 75);
    [Iq, mixing, unmixing, dum, index2centrotypes] = icassoResult(sR,cfg.numcomponent);
    
    % this step is done, because in icassoResult mixing is determined to be
    % pinv(unmixing), which yields strange results. Better take it from the
    % individual iterations. NOTE: as a consequence unmixing*mixing is not
    % necessarily identity anymore !!!
    for k = 1:size(mixing,2)
      ix = sR.index(index2centrotypes(k),:);
      mixing(:,k) = sR.A{ix(1)}(:,ix(2));
    end
    
    %[Iq, mixing, unmixing, dat] = icassoShow(sR, 'estimate', 'off', 'L', cfg.numcomponent);
    
    % sort the output according to Iq
    [srt, ix] = sort(-Iq); % account for NaNs
    mixing    = mixing(:, ix);
    unmixing  = unmixing(ix, :);
    
    cfg.icasso.Iq = Iq(ix);
    cfg.icasso.sR = rmfield(sR, 'signal'); % keep the rest of the information
    
  case 'fastica'
    % check whether the required low-level toolboxes are installed
    ft_hastoolbox('fastica', 1);       % see http://www.cis.hut.fi/projects/ica/fastica
    
    if ~defaultNumCompsUsed &&...
        (~isfield(cfg, 'fastica') || ~isfield(cfg.fastica, 'numOfIC'))
      % user has specified cfg.numcomponent and not specified
      % cfg.fastica.numOfIC, so copy cfg.numcomponent over
      cfg.fastica.numOfIC = cfg.numcomponent;
    elseif ~defaultNumCompsUsed &&...
        isfield(cfg, 'fastica') && isfield(cfg.fastica, 'numOfIC')
      % user specified both cfg.numcomponent and cfg.fastica.numOfIC,
      % unsure which one to use
      ft_error('you can specify either cfg.fastica.numOfIC or cfg.numcomponent (they will have the same effect), but not both');
    end
    
    try
      % construct key-value pairs for the optional arguments
      optarg = ft_cfg2keyval(cfg.fastica);
      [mixing, unmixing] = fastica(dat, optarg{:});
    catch
      % the "catch me" syntax is broken on MATLAB74, this fixes it
      me = lasterror;
      % give a hopefully instructive error message
      ft_info(['If you get an out-of-memory in fastica here, and you use fastica 2.5, change fastica.m, line 482: \n' ...
        'from\n' ...
        '  if ~isempty(W)                  %% ORIGINAL VERSION\n' ...
        'to\n' ...
        '  if ~isempty(W) && nargout ~= 2  %% if nargout == 2, we return [A, W], and NOT ICASIG\n']);
      % forward original error
      rethrow(me);
    end
    
  case 'runica'
    % check whether the required low-level toolboxes are installed
    % see http://www.sccn.ucsd.edu/eeglab
    ft_hastoolbox('eeglab', 1);
    
    if ~defaultNumCompsUsed &&...
        (~isfield(cfg, 'runica') || ~isfield(cfg.runica, 'pca'))
      % user has specified cfg.numcomponent and not specified
      % cfg.runica.pca, so copy cfg.numcomponent over
      cfg.runica.pca = cfg.numcomponent;
    elseif ~defaultNumCompsUsed &&...
        isfield(cfg, 'runica') && isfield(cfg.runica, 'pca')
      % user specified both cfg.numcomponent and cfg.runica.pca,
      % unsure which one to use
      ft_error('you can specify either cfg.runica.pca or cfg.numcomponent (they will have the same effect), but not both');
    end
    
    % construct key-value pairs for the optional arguments
    optarg = [ft_cfg2keyval(cfg.runica) {'reset_randomseed' 0}]; % let FieldTrip deal with the random seed handling
    [weights, sphere] = runica(dat, optarg{:});
    
    % scale the sphering matrix to unit norm
    if strcmp(cfg.normalisesphere, 'yes')
      sphere = sphere./norm(sphere);
    end
    
    unmixing = weights*sphere;
    mixing = [];
    
  case 'binica'
    % check whether the required low-level toolboxes are installed
    % see http://www.sccn.ucsd.edu/eeglab
    ft_hastoolbox('eeglab', 1);
    
    if ~defaultNumCompsUsed &&...
        (~isfield(cfg, 'binica') || ~isfield(cfg.binica, 'pca'))
      % user has specified cfg.numcomponent and not specified
      % cfg.binica.pca, so copy cfg.numcomponent over
      cfg.binica.pca = cfg.numcomponent;
    elseif ~defaultNumCompsUsed &&...
        isfield(cfg, 'binica') && isfield(cfg.binica, 'pca')
      % user specified both cfg.numcomponent and cfg.binica.pca,
      % unsure which one to use
      ft_error('you can specify either cfg.binica.pca or cfg.numcomponent (they will have the same effect), but not both');
    end
    
    % construct key-value pairs for the optional arguments
    optarg = ft_cfg2keyval(cfg.binica);
    [weights, sphere] = binica(dat, optarg{:});
    
    % scale the sphering matrix to unit norm
    if strcmp(cfg.normalisesphere, 'yes')
      sphere = sphere./norm(sphere);
    end
    
    unmixing = weights*sphere;
    mixing = [];
    
  case 'jader'
    % check whether the required low-level toolboxes are installed
    % see http://www.sccn.ucsd.edu/eeglab
    ft_hastoolbox('eeglab', 1);
    
    unmixing = jader(dat, cfg.numcomponent);
    mixing = [];
    
  case 'varimax'
    % check whether the required low-level toolboxes are installed
    % see http://www.sccn.ucsd.edu/eeglab
    ft_hastoolbox('eeglab', 1);
    
    unmixing = varimax(dat);
    mixing = [];
    
  case 'cca'
    % check whether the required low-level toolboxes are installed
    % see http://www.sccn.ucsd.edu/eeglab
    ft_hastoolbox('cca', 1);
    
    [y, w] = ccabss(dat);
    unmixing = w';
    mixing = [];
    
  case 'pca'
    % compute data cross-covariance matrix
    C = (dat*dat')./(size(dat,2)-1);
    
    % eigenvalue decomposition (EVD)
    [E,D] = eig(C);
    
    % sort eigenvectors in descending order of eigenvalues
    d = cat(2,(1:1:Nchans)',diag(D));
    d = sortrows(d, -2);
    
    % return the desired number of principal components
    unmixing = E(:,d(1:cfg.numcomponent,1))';
    mixing = [];
    
    clear C D E d
    
  case 'kpca'
    
    % linear kernel (same as normal covariance)
    %kern = @(X,y) (sum(bsxfun(@times, X, y),2));
    
    % polynomial kernel degree 2
    %kern = @(X,y) (sum(bsxfun(@times, X, y),2).^2);
    
    % RBF kernel
    kern = @(X,y) (exp(-0.5* sqrt(sum(bsxfun(@minus, X, y).^2, 2))));
    
    % compute kernel matrix
    C = zeros(Nchans,Nchans);
    ft_progress('init', cfg.feedback, 'computing kernel matrix...');
    for k = 1:Nchans
      ft_progress(k/Nchans, 'computing kernel matrix %d from %d', k, Nchans);
      C(k,:) = kern(dat, dat(k,:));
    end
    ft_progress('close');
    
    % eigenvalue decomposition (EVD)
    [E,D] = eig(C);
    
    % sort eigenvectors in descending order of eigenvalues
    d = cat(2,(1:1:Nchans)',diag(D));
    d = sortrows(d, -2);
    
    % return the desired number of principal components
    unmixing = E(:,d(1:cfg.numcomponent,1))';
    mixing = [];
    
    clear C D E d
    
  case 'svd'
    % it is more memory efficient to use the (non-scaled) covariance
    if cfg.numcomponent<Nchans
      % compute only the first components
      [u, s, v] = svds(dat*dat', cfg.numcomponent);
    else
      % compute all components
      [u, s, v] = svd(dat*dat', 0);
    end
    clear s v % not needed
    
    unmixing = u';
    mixing = [];
    
  case 'dss'
    % check whether the required low-level toolboxes are installed
    % see http://www.cis.hut.fi/projects/dss
    ft_hastoolbox('dss', 1);
    
    params         = struct(cfg.dss);
    params.denf.h  = str2func(cfg.dss.denf.function);
    params.preprocf.h = str2func(cfg.dss.preprocf.function);
    if ~ischar(cfg.numcomponent)
      params.sdim = cfg.numcomponent;
    end
    if isfield(cfg.dss, 'wdim') && ~isempty(cfg.dss.wdim)
      params.wdim = cfg.dss.wdim;
    end
    if isfield(cfg.dss, 'V') && ~isempty(cfg.dss.V)
      params.Y = params.V*dat;
    end
    
    % create the state
    state   = dss_create_state(dat, params);
    if isfield(cfg.dss, 'V') && ~isempty(cfg.dss.V)
      state.V = cfg.dss.V;
    end
    if isfield(cfg.dss, 'dV') && ~isempty(cfg.dss.dV)
      state.dV = cfg.dss.dV;
    end
    
    % increase the amount of information that is displayed on screen
    % state.verbose = 3;
    % start the decomposition
    % state   = dss(state);  % this is for the DSS toolbox version 0.6 beta
    state   = denss(state);  % this is for the DSS toolbox version 1.0
    
    mixing   = state.A;
    unmixing = state.B;
  
    % remember the updated configuration details
    cfg.dss.denf      = state.denf;
    cfg.dss.orthof    = state.orthof;
    cfg.dss.preprocf  = state.preprocf;
    cfg.dss.stopf     = state.stopf;
    cfg.dss.W         = state.W;
    cfg.dss.V         = state.V;
    cfg.dss.dV        = state.dV;
    if isfield(state, 'D'), cfg.dss.D = state.D(1:min([state.sdim size(state.dV)])); end
    cfg.numcomponent  = min([state.sdim size(state.dV)]);
    
  case 'sobi'
    % check whether the required low-level toolboxes are installed
    % see http://www.sccn.ucsd.edu/eeglab
    ft_hastoolbox('eeglab', 1);
    
    % check for additional options, see SOBI for details
    if ~isfield(cfg, 'sobi')
      mixing = sobi(dat, cfg.numcomponent);
    elseif isfield(cfg.sobi, 'n_sources') && isfield(cfg.sobi, 'p_correlations')
      mixing = sobi(dat, cfg.sobi.n_sources, cfg.sobi.p_correlations);
    elseif isfield(cfg.sobi, 'n_sources')
      mixing = sobi(dat,cfg.sobi.n_sources);
    else
      ft_error('unknown options for SOBI component analysis');
    end
    
    unmixing = [];
    
  case 'predetermined unmixing matrix'
    % check which labels from the cfg are identical to those of the data
    % this gives us the rows of cfg.topo (the channels) and of
    % data.trial (also channels) that we are going to use later
    [datsel, chansel] = match_str(data.label, cfg.topolabel);
    
    % ensure 1:1 corresponcence between cfg.topolabel & data.label
    % otherwise we cannot compute the components (if source channels are
    % missing) or will have a problem when projecting it back (because we
    % dont have a marker to say that there are channels in data.label
    % which we did not use and thus can't recover from source-space)
    
    if length(cfg.topolabel)<length(chansel)
      ft_error('cfg.topolabels do not uniquely correspond to data.label, please check')
    end
    if length(data.label)<length(datsel)
      ft_error('cfg.topolabels do not uniquely correspond to data.label, please check')
    end
    
    % reorder the mixing matrix so that the channel order matches the order in the data
    cfg.unmixing  = cfg.unmixing(:,chansel);
    cfg.topolabel = cfg.topolabel(chansel);
    
    unmixing = cfg.unmixing;
    mixing   = [];
    
  case 'white'
    % compute the covariance matrix and an unmixing matrix that makes the data white
    c = dat*dat';
    c = c./(size(dat,2)-1);
    [u, s] = svd(c);
    % split the singular values into half
    for i=1:size(s)
      if (s(i,i)/s(1,1))>(100*eps)
        s(i,i) = 1./sqrt(s(i,i));
      else
        s(i,i) = 0;
      end
    end
    unmixing = s * u';
    mixing   = [];
    
  case 'csp'
    C1 = cov(dat1');
    C2 = cov(dat2');
    unmixing = csp(C1, C2, cfg.csp.numfilters);
    mixing   = [];  % will be computed below
    
  case 'bsscca'
    % this method relies on time shifting of the original data, in much the
    % same way as ft_denoise_tsr. as such it is more natural to represent
    % the data in the cell-array, because the trial-boundaries are clear.
    % if represented in a concatenated array one has to keep track of the
    % trial boundaries
    
    optarg          = ft_cfg2keyval(cfg.bsscca);
    optarg          = cat(2,optarg, {'time', data.time});
    [unmixing, mixing, rho, compdata, time] = bsscca(dat, optarg{:});
    data.trial = mixing*compdata;
    data.time  = time;
    
    if size(mixing,1)>numel(data.label)
      for m = 1:(size(mixing,1)-numel(data.label))
        data.label{end+1} = sprintf('refchan%03d',m);
      end
    end
    
    % remember the canonical correlations
    cfg.bsscca.rho = rho;
    
  case 'parafac'
    ft_error('parafac is not supported anymore in ft_componentanalysis');
    
  otherwise
    ft_error('unknown method for component analysis');
end % switch method

% make sure we have both mixing and unmixing matrices
% if not, compute (pseudo-)inverse to go from one to the other
if isempty(unmixing) && ~isempty(mixing)
  if (size(mixing,1)==size(mixing,2))
    unmixing = inv(mixing);
  else
    unmixing = pinv(mixing);
  end
elseif isempty(mixing) && ~isempty(unmixing)
  if (size(unmixing,1)==size(unmixing,2)) && rank(unmixing)==size(unmixing,1)
    mixing = inv(unmixing);
  else
    mixing = pinv(unmixing);
  end
elseif isempty(mixing) && isempty(unmixing)
  % this sanity check is needed to catch convergence problems in fastica
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1519
  ft_error('the component unmixing failed');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results and construct data structure

comp = [];
if isfield(data, 'fsample'), comp.fsample = data.fsample; end
if isfield(data, 'time'),    comp.time    = data.time;    end

% make sure we don't return more components than were requested
% (some methods respect the maxcomponent parameters, others just always
% return a fixed (i.e., numchans) number of components)
if size(unmixing,1) > cfg.numcomponent
  unmixing(cfg.numcomponent+1:end,:) = [];
end
if size(mixing,2) > cfg.numcomponent
  mixing(:,cfg.numcomponent+1:end) = [];
end

% compute the activations in each trial
if strcmp(cfg.doscale, 'yes')
  for trial=1:Ntrials
    comp.trial{trial} = scale * unmixing * data.trial{trial};
  end
else
  for trial=1:Ntrials
    comp.trial{trial} = unmixing * data.trial{trial};
  end
end

% store mixing/unmixing matrices in structure
comp.topo = mixing;
comp.unmixing = unmixing;

% get the labels
if strcmp(cfg.method, 'predetermined unmixing matrix')
  prefix = 'component';
else
  prefix = cfg.method;
end

for k = 1:size(comp.topo,2)
  comp.label{k,1} = sprintf('%s%03d', prefix, k);
end
comp.topolabel = data.label(:);

if isfield(data, 'grad')
  sensfield = 'grad';
elseif isfield(data, 'elec')
  sensfield = 'elec';
elseif isfield(data, 'opto')
  sensfield = 'opto';
else
  sensfield = [];
end

% apply the linear projection also to the sensor description
if ~isempty(sensfield)
  if  strcmp(cfg.updatesens, 'yes')
    ft_info('also applying the unmixing matrix to the %s structure\n', sensfield);
    % construct a montage and apply it to the sensor description
    montage          = [];
    montage.labelold = data.label;
    montage.labelnew = comp.label;
    montage.tra      = unmixing;
    comp.(sensfield) = ft_apply_montage(data.(sensfield), montage, 'balancename', 'comp', 'keepunused', 'yes');
    
    % The output sensor array cannot simply be interpreted as the input
    % sensor array, hence the type should be removed to allow autodetection
    % See also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1806
    if isfield(comp.(sensfield), 'type')
      comp.(sensfield) = rmfield(comp.(sensfield), 'type');
    end
  else
    ft_info('not applying the unmixing matrix to the %s structure\n', sensfield);
    % simply copy it over
    comp.(sensfield) = data.(sensfield);
  end
end % if sensfield

% copy the sampleinfo into the output
if isfield(data, 'sampleinfo')
  comp.sampleinfo = data.sampleinfo;
end

% copy the trialinfo into the output
if isfield(data, 'trialinfo')
  comp.trialinfo = data.trialinfo;
end

% convert back to input type if necessary
if istimelock
  comp = ft_checkdata(comp, 'datatype', 'timelock+comp');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble randomseed
ft_postamble previous   data
ft_postamble provenance comp
ft_postamble history    comp
ft_postamble savevar    comp
