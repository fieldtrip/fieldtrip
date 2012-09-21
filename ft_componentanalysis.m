function [comp] = ft_componentanalysis(cfg, data)

% FT_COMPONENTANALYSIS performs independent component analysis or other
% spatio-temporal decompositions of EEG or MEG data. This function computes
% the topography and timecourses of the components. The output of this
% function can be further analyzed with FT_TIMELOCKANALYSIS or
% FT_FREQNANALYSIS.
%
% Use as
%   [comp] = ft_componentanalysis(cfg, data)
%
% where the data comes from FT_PREPROCESSING and the configuration
% structure can contain
%   cfg.method       = 'runica', 'fastica', 'binica', 'pca', 'svd', 'jader', 'varimax', 'dss', 'cca', 'sobi', 'white' or 'csp' (default = 'runica')
%   cfg.channel      = cell-array with channel selection (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.trials       = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.numcomponent = 'all' or number (default = 'all')
%   cfg.demean       = 'no' or 'yes' (default = 'yes')
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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
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
ft_preamble randomseed
ft_preamble trackconfig
ft_preamble loadvar data

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', {'detrend'});
cfg = ft_checkconfig(cfg, 'renamed',   {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamedval',   {'method', 'predetermined mixing matrix', 'predetermined unmixing matrix'});
cfg = ft_checkconfig(cfg, 'deprecated', {'topo'});

% set the defaults
cfg.method          = ft_getopt(cfg, 'method',       'runica');
cfg.demean          = ft_getopt(cfg, 'demean',       'yes');
cfg.trials          = ft_getopt(cfg, 'trials',       'all');
cfg.channel         = ft_getopt(cfg, 'channel',      'all');
cfg.numcomponent    = ft_getopt(cfg, 'numcomponent', 'all');
cfg.inputfile       = ft_getopt(cfg, 'inputfile',    []);
cfg.outputfile      = ft_getopt(cfg, 'outputfile',   []);
cfg.normalisesphere = ft_getopt(cfg, 'normalisesphere', 'yes');

% select channels, has to be done prior to handling of previous (un)mixing matrix
cfg.channel = ft_channelselection(cfg.channel, data.label);

if isfield(cfg, 'topo') && isfield(cfg, 'topolabel')
  warning(['Specifying cfg.topo (= mixing matrix) to determine component '...
    'timecourses in specified data is deprecated; please specify an '...
    'unmixing matrix instead with cfg.unmixing. '...
    'Using cfg.unmixing=pinv(cfg.topo) for now to reproduce old behaviour.']);
  
  cfg.unmixing = pinv(cfg.topo);
  cfg = rmfield(cfg,'topo');
end

if isfield(cfg, 'unmixing') && isfield(cfg, 'topolabel')
  % use the previously determined unmixing matrix on this dataset
  
  % test whether all required channels are present in the data
  [datsel, toposel] = match_str(cfg.channel, cfg.topolabel);
  if length(toposel)~=length(cfg.topolabel)
    error('not all channels that are required for the unmixing are present in the data');
  end
  
  % ensure that all data channels not used in the unmixing should be removed from the channel selection
  tmpchan = match_str(cfg.channel, cfg.topolabel);
  cfg.channel = cfg.channel(tmpchan);
  
  % remove all cfg settings  that do not apply
  tmpcfg              = [];
  tmpcfg.demean       = cfg.demean;
  tmpcfg.trials       = cfg.trials;
  tmpcfg.unmixing     = cfg.unmixing;    % the NxM unmixing matrix (M channels, N components)
  tmpcfg.topolabel    = cfg.topolabel;   % the Mx1 labels of the data that was used in determining the mixing matrix
  tmpcfg.channel      = cfg.channel;     % the Mx1 labels of the data that is presented now to this function
  tmpcfg.numcomponent = 'all';
  tmpcfg.method       = 'predetermined unmixing matrix';
  tmpcfg.outputfile   = cfg.outputfile;
  cfg                 = tmpcfg;
end

% add the options for the specified methods to the configuration, only if needed
switch cfg.method
case 'icasso'
  cfg.icasso        = ft_getopt(cfg,        'icasso', []);
  cfg.icasso.mode   = ft_getopt(cfg.icasso, 'mode',   'both');
  cfg.icasso.Niter  = ft_getopt(cfg.icasso, 'Niter',  15);
  cfg.icasso.method = 'fastica'; % as of yet only this one is supported
  
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
case 'csp'
  % additional options, see CSP for details
  cfg.csp = ft_getopt(cfg, 'csp', []);
  cfg.csp.numfilters = ft_getopt(cfg.csp, 'numfilters', 6);
  cfg.csp.classlabels = ft_getopt(cfg.csp, 'classlabels');
case 'bsscca'
  % additional options, see BSSCCA for details
  cfg.bsscca       = ft_getopt(cfg,        'bsscca', []);
  cfg.bsscca.delay = ft_getopt(cfg.bsscca, 'delay', 1);
otherwise
  % do nothing
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = ft_selectdata(data, 'rpt', cfg.trials);
end
Ntrials  = length(data.trial);

% select channels of interest
chansel = match_str(data.label, cfg.channel);
fprintf('selecting %d channels\n', length(chansel));
for trial=1:Ntrials
  data.trial{trial} = data.trial{trial}(chansel,:);
end
data.label = data.label(chansel);
Nchans     = length(chansel);
if Nchans==0
  error('no channels were selected');
end

% default is to compute just as many components as there are channels in the data
if strcmp(cfg.numcomponent, 'all')
  cfg.numcomponent = length(data.label);
end

% determine the size of each trial, they can be variable length
Nsamples = zeros(1,Ntrials);
for trial=1:Ntrials
  Nsamples(trial) = size(data.trial{trial},2);
end

if strcmp(cfg.demean, 'yes')
  % optionally perform baseline correction on each trial
  fprintf('baseline correcting data \n');
  for trial=1:Ntrials
    data.trial{trial} = ft_preproc_baselinecorrect(data.trial{trial});
  end
end

% determine the scaling of the data, scale it to approximately unity
% this will improve the performance of some methods, esp. fastica
scale = norm((data.trial{1}*data.trial{1}')./size(data.trial{1},2));
scale = sqrt(scale);
if scale ~= 0
  fprintf('scaling data with 1 over %f\n', scale);
  for trial=1:Ntrials
    data.trial{trial} = data.trial{trial} ./ scale;
  end
else
  fprintf('no scaling applied, since factor is 0\n');
end

if strcmp(cfg.method, 'sobi')
  
  % concatenate all the data into a 3D matrix respectively 2D (sobi)
  fprintf('concatenating data');
  Nsamples = Nsamples(1);
  dat = zeros(Ntrials, Nchans, Nsamples);
  % all trials should have an equal number of samples
  % and it is assumed that the time axes of all trials are aligned
  for trial=1:Ntrials
    fprintf('.');
    dat(trial,:,:) = data.trial{trial};
  end
  fprintf('\n');
  fprintf('concatenated data matrix size %dx%dx%d\n', size(dat,1), size(dat,2), size(dat,3));
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
    error('CSP requires class labels!');
  end
  if length(sel1)+length(sel2)~=length(cfg.csp.classlabels)
    warning('not all trials belong to class 1 or 2');
  end
  dat1 = cat(2, data.trial{sel1});
  dat2 = cat(2, data.trial{sel2});
  fprintf('concatenated data matrix size for class 1 is %dx%d\n', size(dat1,1), size(dat1,2));
  fprintf('concatenated data matrix size for class 2 is %dx%d\n', size(dat2,1), size(dat2,2));
  
elseif (~strcmp(cfg.method, 'predetermined unmixing matrix') && ~strcmp(cfg.method, 'bsscca')) && ~strcmp(cfg.method, 'dss2')
  
  % concatenate all the data into a 2D matrix unless we already have an
  % unmixing matrix
  fprintf('concatenating data');
  
  dat = zeros(Nchans, sum(Nsamples));
  for trial=1:Ntrials
    fprintf('.');
    begsample = sum(Nsamples(1:(trial-1))) + 1;
    endsample = sum(Nsamples(1:trial));
    dat(:,begsample:endsample) = data.trial{trial};
  end
  fprintf('\n');
  fprintf('concatenated data matrix size %dx%d\n', size(dat,1), size(dat,2));
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform the component analysis
fprintf('starting decomposition using %s\n', cfg.method);
switch cfg.method
 
  case 'icasso'
    % check whether the required low-level toolboxes are installed
    ft_hastoolbox('icasso',  1);
    
    if strcmp(cfg.icasso.method, 'fastica')
      ft_hastoolbox('fastica', 1); 
      cfg.fastica.numOfIC = cfg.numcomponent;  
    else
      error('only ''fastica'' is supported as method for icasso');
    end

    optarg = ft_cfg2keyval(cfg.(cfg.icasso.method));
    sR     = icassoEst(cfg.icasso.mode, dat, cfg.icasso.Niter, optarg{:});
    sR     = icassoExp(sR);
    [Iq, mixing, unmixing] = icassoShow(sR);    

    cfg.icasso.Iq = Iq;
    cfg.icasso.sR = rmfield(sR, 'signal');

  case 'fastica'
    % check whether the required low-level toolboxes are installed
    ft_hastoolbox('fastica', 1);       % see http://www.cis.hut.fi/projects/ica/fastica
    
    % set the number of components to be estimated
    cfg.fastica.numOfIC = cfg.numcomponent;
    
    try
      % construct key-value pairs for the optional arguments
      optarg = ft_cfg2keyval(cfg.fastica);
      [mixing, unmixing] = fastica(dat, optarg{:});
    catch
      % the "catch me" syntax is broken on MATLAB74, this fixes it
      me = lasterror;
      % give a hopefully instructive error message
      fprintf(['If you get an out-of-memory in fastica here, and you use fastica 2.5, change fastica.m, line 482: \n' ...
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
    
    % construct key-value pairs for the optional arguments
    optarg = ft_cfg2keyval(cfg.runica);
    [weights, sphere] = runica(dat, optarg{:});
    
    % scale the sphering matrix to unit norm
    if strcmp(cfg.normalisesphere, 'yes'),
      sphere = sphere./norm(sphere);
    end
    
    unmixing = weights*sphere;
    mixing = [];
    
  case 'binica'
    % check whether the required low-level toolboxes are installed
    % see http://www.sccn.ucsd.edu/eeglab
    ft_hastoolbox('eeglab', 1);
    
    % construct key-value pairs for the optional arguments
    optarg = ft_cfg2keyval(cfg.binica);
    [weights, sphere] = binica(dat, optarg{:});
    
    % scale the sphering matrix to unit norm
    if strcmp(cfg.normalisesphere, 'yes'),
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
    d = sortrows(d,[-2]);
    
    % return the desired number of principal components
    unmixing = E(:,d(1:cfg.numcomponent,1))';
    mixing = [];
    
    clear C D E d
    
  case 'svd'
    if cfg.numcomponent<Nchans
      % compute only the first components
      [u, s, v] = svds(dat, cfg.numcomponent);
    else
      % compute all components
      [u, s, v] = svd(dat, 0);
    end
    
    unmixing = u';
    mixing = [];
    
  case 'dss'
    % check whether the required low-level toolboxes are installed
    % see http://www.cis.hut.fi/projects/dss
    ft_hastoolbox('dss', 1);
    
    params         = struct(cfg.dss);
    params.denf.h  = str2func(cfg.dss.denf.function);
    if ~ischar(cfg.numcomponent)
      params.sdim = cfg.numcomponent;
    end
    % create the state
    state   = dss_create_state(dat, params);
    % increase the amount of information that is displayed on screen
    %state.verbose = 3;
    % start the decomposition
    % state   = dss(state);  % this is for the DSS toolbox version 0.6 beta
    state   = denss(state);  % this is for the DSS toolbox version 1.0
    %weights = state.W;
    %sphere  = state.V;
    
    mixing   = state.A;
    unmixing = state.B;
    
    % remember the updated configuration details
    cfg.dss.denf      = state.denf;
    cfg.dss.orthof    = state.orthof;
    cfg.dss.preprocf  = state.preprocf;
    cfg.dss.stopf     = state.stopf;
    cfg.dss.W         = state.W;
    cfg.dss.V         = state.V;
    cfg.numcomponent  = state.sdim;
    
  case 'dss2'
    % check whether the required low-level toolboxes are installed
    % see http://www.cis.hut.fi/projects/dss
    %ft_hastoolbox('dss2', 1);
    
    params         = struct(cfg.dss);
    params.denf.h  = str2func(cfg.dss.denf.function);
    if ~ischar(cfg.numcomponent)
      params.sdim = cfg.numcomponent;
    end
    % create the state
    state   = dss_create_state(data.trial, params);
    % increase the amount of information that is displayed on screen
    % state.verbose = 3;
    % start the decomposition
    % state   = dss(state);  % this is for the DSS toolbox version 0.6 beta
    state   = denss(state);  % this is for the DSS toolbox version 1.0
    %weights = state.W;
    %sphere  = state.V;
    
    mixing   = state.A;
    unmixing = state.B;
    
    % remember the updated configuration details
    cfg.dss.denf      = state.denf;
    cfg.dss.orthof    = state.orthof;
    cfg.dss.preprocf  = state.preprocf;
    cfg.dss.stopf     = state.stopf;
    cfg.dss.W         = state.W;
    cfg.dss.V         = state.V;
    cfg.numcomponent  = state.sdim;
    
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
      error('unknown options for SOBI component analysis');
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
      error('COMPONENTANALYSIS:LABELMISSMATCH:topolabel', 'cfg.topolabels do not uniquely correspond to data.label, please check')
    end
    if length(data.label)<length(datsel)
      error('COMPONENTANALYSIS:LABELMISSMATCH:topolabel', 'cfg.topolabels do not uniquely correspond to data.label, please check')
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
    
    
    [unmixing, rho] = bsscca(data.trial,cfg.bsscca.delay);
    mixing          = [];
    %unmixing        = diag(rho);
    
  case 'parafac'
    error('parafac is not supported anymore in ft_componentanalysis');
    
  otherwise
    error('unknown method for component analysis');
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
  % see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1519
  error('the component unmixing failed');
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
for trial=1:Ntrials
  comp.trial{trial} = scale * unmixing * data.trial{trial};
end

% store mixing/unmixing matrices in structure
comp.topo = mixing;
comp.unmixing = unmixing;

% get the labels
if strcmp(cfg.method, 'predetermined unmixing matrix'),
  prefix = 'component';
else
  prefix = cfg.method;
end

for k = 1:size(comp.topo,2)
  comp.label{k,1} = sprintf('%s%03d', prefix, k);
end
comp.topolabel = data.label(:);

% apply the montage also to the elec/grad, if present
if isfield(data, 'grad') || (isfield(data, 'elec') && isfield(data.elec, 'tra'))
  fprintf('applying the mixing matrix to the sensor description\n');
  if isfield(data, 'grad')
    sensfield = 'grad';
  else
    sensfield = 'elec';
  end
  montage          = [];
  montage.labelorg = data.label;
  montage.labelnew = comp.label;
  montage.tra      = unmixing;
  comp.(sensfield) = ft_apply_montage(data.(sensfield), montage, 'balancename', 'comp', 'keepunused', 'yes');
end

% copy the sampleinfo into the output
if isfield(data, 'sampleinfo')
  comp.sampleinfo = data.sampleinfo;
end

% copy the trialinfo into the output
if isfield(data, 'trialinfo')
  comp.trialinfo = data.trialinfo;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble randomseed
ft_postamble previous data
ft_postamble history comp
ft_postamble savevar comp
