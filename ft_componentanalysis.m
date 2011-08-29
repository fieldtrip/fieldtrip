function [comp] = ft_componentanalysis(cfg, data)

% FT_COMPONENTANALYSIS principal or independent component analysis
% computes the topography and timecourses of the ICA/PCA components
% in the EEG/MEG data.
%
% Use as
%   [comp] = ft_componentanalysis(cfg, data)
%
% where the data comes from FT_PREPROCESSING or FT_TIMELOCKANALYSIS and the
% configuration structure can contain
%   cfg.method       = 'runica', 'fastica', 'binica', 'pca', 'jader', 'varimax', 'dss', 'cca', 'sobi' (default = 'runica')
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
% Instead of specifying a component analysis method, you can also specify
% a previously computed mixing matrix, which will be used to estimate the
% component timecourses in this data. This requires
%   cfg.topo         = NxN matrix with a component topography in each column
%   cfg.topolabel    = Nx1 cell-array with the channel labels
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
% See also FASTICA, RUNICA, BINICA, SVD, JADER, VARIMAX, DSS, CCA, SOBI

% NOTE parafac is also implemented, but that does not fit into the
% structure of 2D decompositions very well. Probably I should implement it
% in a separate function for N-D decompositions

% Copyright (C) 2003-2011, Robert Oostenveld
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

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();

% set a timer to determine how long this function takes
stopwatch=tic;

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'forbidden', {'detrend'});
cfg = ft_checkconfig(cfg, 'renamed',   {'blc', 'demean'});

% set the defaults
cfg.method       = ft_getopt(cfg, 'method',       'runica');
cfg.demean       = ft_getopt(cfg, 'demean',       'yes');
cfg.trials       = ft_getopt(cfg, 'trials',       'all');
cfg.channel      = ft_getopt(cfg, 'channel',      'all');
cfg.numcomponent = ft_getopt(cfg, 'numcomponent', 'all');
cfg.inputfile    = ft_getopt(cfg, 'inputfile',    []);
cfg.outputfile   = ft_getopt(cfg, 'outputfile',   []);
cfg.normalisesphere = ft_getopt(cfg, 'normalisesphere', 'yes');

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

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% select channels, has to be done prior to handling of previous (un)mixing matrix
cfg.channel = ft_channelselection(cfg.channel, data.label);

if isfield(cfg, 'topo') && isfield(cfg, 'topolabel')
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
  tmpcfg.topo         = cfg.topo;        % the MxN mixing matrix (M channels, N components)
  tmpcfg.topolabel    = cfg.topolabel;   % the Mx1 labels of the data that was used in determining the mixing matrix
  tmpcfg.channel      = cfg.channel;     % the Mx1 labels of the data that is presented now to this function
  tmpcfg.numcomponent = 'all';
  tmpcfg.method       = 'predetermined mixing matrix';
  tmpcfg.outputfile   = cfg.outputfile;
  cfg                 = tmpcfg;
end

% additional options, see FASTICA for details
cfg.fastica = ft_getopt(cfg, 'fastica', []);

% additional options, see RUNICA for details
cfg.runica       = ft_getopt(cfg,        'runica',  []);
cfg.runica.lrate = ft_getopt(cfg.runica, 'lrate',   0.001);

% additional options, see BINICA for details
cfg.binica       = ft_getopt(cfg,        'binica',  []);
cfg.binica.lrate = ft_getopt(cfg.binica, 'lrate',   0.001);

% additional options, see DSS for details
cfg.dss               = ft_getopt(cfg,          'dss',      []);
cfg.dss.denf          = ft_getopt(cfg.dss,      'denf',     []);
cfg.dss.denf.function = ft_getopt(cfg.dss.denf, 'function', 'denoise_fica_tanh');
cfg.dss.denf.params   = ft_getopt(cfg.dss.denf, 'params',   []);

% check whether the required low-level toolboxes are installed
switch cfg.method
  case 'fastica'
    ft_hastoolbox('fastica', 1);       % see http://www.cis.hut.fi/projects/ica/fastica
  case {'runica', 'jader', 'varimax', 'binica', 'sobi'}
    ft_hastoolbox('eeglab', 1);        % see http://www.sccn.ucsd.edu/eeglab
  case 'parafac'
    ft_hastoolbox('nway', 1);          % see http://www.models.kvl.dk/source/nwaytoolbox
  case 'dss'
    ft_hastoolbox('dss', 1);           % see http://www.cis.hut.fi/projects/dss
end % cfg.method

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

if strcmp(cfg.method, 'predetermined mixing matrix')
  % the single trial data does not have to be concatenated
elseif strcmp(cfg.method, 'parafac') || strcmp(cfg.method, 'sobi')
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
  if strcmp(cfg.method, 'sobi')
    % sobi wants Nchans, Nsamples, Ntrials matrix and for Ntrials = 1 remove
    % trial dimension
    if Ntrials == 1
        dummy = 0;
        [dat, dummy] = shiftdim(dat);
    else
        dat = shiftdim(dat,1); 
    end
  end 
else
  % concatenate all the data into a 2D matrix
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
  
  case 'fastica'
    
    try
      % construct key-value pairs for the optional arguments
      optarg = ft_cfg2keyval(cfg.fastica);
      [A, W] = fastica(dat, optarg{:});
      weights = W;
      sphere = eye(size(W,2));
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
    % construct key-value pairs for the optional arguments
    optarg = ft_cfg2keyval(cfg.runica);
    [weights, sphere] = runica(dat, optarg{:});

    % scale the sphering matrix to unit norm
    if strcmp(cfg.normalisesphere, 'yes'),
      sphere = sphere./norm(sphere);
    end
    
  case 'binica'
    % construct key-value pairs for the optional arguments
    optarg = ft_cfg2keyval(cfg.binica);
    [weights, sphere] = binica(dat, optarg{:});
    
    % scale the sphering matrix to unit norm
    if strcmp(cfg.normalisesphere, 'yes'),
      sphere = sphere./norm(sphere);
    end
    
  case 'jader'
    weights = jader(dat);
    sphere  = eye(size(weights, 2));
    
  case 'varimax'
    weights = varimax(dat);
    sphere  = eye(size(weights, 2));
    
  case 'cca'
    [y, w] = ccabss(dat);
    weights = w';
    sphere  = eye(size(weights, 2));
    
  case 'pca'
    % compute data cross-covariance matrix
    C = (dat*dat')./(size(dat,2)-1);
    % eigenvalue decomposition (EVD)
    [E,D] = eig(C);
    % sort eigenvectors in descending order of eigenvalues
    d = cat(2,(1:1:Nchans)',diag(D));
    d = sortrows(d,[-2]);
    % return the desired number of principal components
    weights = E(:,d(1:cfg.numcomponent,1))';
    sphere = eye(size(weights,2));
    clear C D E d
    
  case 'svd'
    if cfg.numcomponent<Nchans
      % compute only the first components
      [u, s, v] = svds(dat, cfg.numcomponent);
      u(Nchans, Nchans) = 0;
    else
      % compute all components
      [u, s, v] = svd(dat, 0);
    end
    weights = u';
    sphere  = eye(size(weights, 2));
    
  case 'parafac'
    f = parafac(dat, cfg.numcomponent);
    
  case 'dss'
    params         = struct(cfg.dss);
    params.denf.h  = str2func(cfg.dss.denf.function);
    if ~ischar(cfg.numcomponent)
      params.sdim = cfg.numcomponent;
    end
    % create the state
    state   = dss_create_state(dat, params);
    % increase the amount of information that is displayed on screen
    state.verbose = 3;
    % start the decomposition
    % state   = dss(state);  % this is for the DSS toolbox version 0.6 beta
    state   = denss(state);  % this is for the DSS toolbox version 1.0
    %weights = state.W;
    %sphere  = state.V;
    if size(state.A,1)==size(state.A,2)
      weights  = inv(state.A);
    else
      weights  = pinv(state.A);
    end
    sphere   = eye(size(weights, 2));
    % remember the updated configuration details
    cfg.dss.denf      = state.denf;
    cfg.numcomponent  = state.sdim;

  case 'sobi'
    % check for additional options, see SOBI for details
    if ~isfield(cfg, 'sobi')
        mixm = sobi(dat);
    elseif isfield(cfg.sobi, 'n_sources') && isfield(cfg.sobi, 'p_correlations')
        mixm = sobi(dat, cfg.sobi.n_sources, cfg.sobi.p_correlations);
    elseif isfield(cfg.sobi, 'n_sources')
        mixm = sobi(dat,cfg.sobi.n_sources);
    else
        error('unknown options for SOBI component analysis');      
    end
    weights = pinv(mixm);
    sphere  = eye(size(weights, 2));

  case 'predetermined mixing matrix'
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
    cfg.topo      = cfg.topo(chansel,:);
    cfg.topolabel = cfg.topolabel(chansel);
    
    % get the number of channels we are using from the data
    Nchan   = size(cfg.topo, 1);
    % get the number of components in which the data was decomposed
    
    Ncomp   = size(cfg.topo, 2);
    cfg.numcomponent = Ncomp;
    
    % initialize sphere and weights
    sphere  = eye(Nchan,Nchan);
    
    % cfg.topo is a Channel x Component matrix (the mixing matrix, A)
    % lets get the unmixing matrix (weights, W)
    
    % now, the weights matrix is simply given by its (pseudo)-inverse
    if (Nchan==Ncomp)
      weights = inv(cfg.topo);
    else
      weights = pinv(cfg.topo);
    end
    
  otherwise
    error('unknown method for component analysis');
end % switch method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
if isfield(data, 'fsample'), comp.fsample = data.fsample; end
if isfield(data, 'time'),    comp.time    = data.time;    end

% compute the activations in each trial
for trial=1:Ntrials
  if strcmp(cfg.method, 'parafac')
    % FIXME, this is not properly supported yet
    comp.trial{trial} = [];
  else
    comp.trial{trial} = weights * sphere * data.trial{trial};
  end
end

% get the mixing matrix
if strcmp(cfg.method, 'parafac')
  comp.topo = f{2};
  comp.f1   = f{1}; %FIXME, this is not properly supported yet
  comp.f2   = f{2};
  comp.f3   = f{3};
elseif size(weights,1)==size(weights,2)
  comp.topo = inv(weights*sphere);
else
  comp.topo = pinv(weights*sphere); %allow fewer sources than sensors
end

% get the labels
if strcmp(cfg.method, 'predetermined mixing matrix'),
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
  montage.tra      = weights * sphere;
  comp.(sensfield) = ft_apply_montage(data.(sensfield), montage, 'balancename', 'comp', 'keepunused', 'yes');
end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id   = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();

% remember the configuration details of the input data
if isfield(data, 'cfg'), cfg.previous = data.cfg; end

% remember the exact configuration details in the output
comp.cfg = cfg;

% copy the sampleinfo into the output
if isfield(data, 'sampleinfo')
  comp.sampleinfo = data.sampleinfo;
end
% copy the trialinfo into the output
if isfield(data, 'trialinfo')
  comp.trialinfo = data.trialinfo;
end

fprintf('total time in componentanalysis %.1f seconds\n', toc(stopwatch));

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'comp', comp); % use the variable name "data" in the output file
end
