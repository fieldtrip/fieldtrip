function dataout = ft_regressconfound(cfg, datain)

% FT_REGRESSCONFOUND estimates the regression weight of a set of confounds
% (GLM) and removes the estimated contribution from the single-trial data
%
% Use as
%  timelock = ft_regressconfound(cfg, timelock)
% or
%  freq = ft_regressconfound(cfg, freq)
%
% where datain comes from FT_TIMELOCKANALYSIS or FT_FREQANALYSIS with
% keeptrials = 'yes' and cfg is a configuratioun structure that should
% contain
%
%   cfg.confound    = matrix, [Ntrials X Nconfounds]
%
% The following configuration options are supported:
%   cfg.reject      = vector, [1 X Nconfounds], listing the confounds that
%                     are to be rejected (default = 'all')
%   cfg.normalize   = string, 'yes' or 'no', normalization to
%                     make the confounds orthogonal (default = 'yes')
%   cfg.statistics  = string, 'yes' or 'no', whether to add the statistics
%                     to the output (default = 'no')
%   cfg.model       = string, 'yes' or 'no', whether to add the model to
%                     the output (default = 'no')
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
% See also FT_REJECTCOMPONENT, FT_REJECTARTIFACT

% Copyrights (C) 2011, Robert Oostenveld, Arjen Stolk, Lennart Verhagen
%
% $Id: ft_rejectconfound.m 4375 2011-10-07 09:16:45Z arjsto $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();
ftFuncMem   = memtic();

% defaults for optional input/ouputfile
cfg.inputfile  = ft_getopt(cfg, 'inputfile',  []);
cfg.outputfile = ft_getopt(cfg, 'outputfile', []);

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    datain = loadvar(cfg.inputfile, 'data');
  end
end

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
datain = ft_checkdata(datain, 'datatype', {'timelock', 'freq'}, 'feedback', 'yes', 'hastrials', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'confound'});

% confound specification
regr = ft_getopt(cfg, 'confound');  % there is no default value
nconf = size(regr,2);
conflist = 1:nconf;
if ~isfield(cfg, 'reject') || strcmp(cfg.reject, 'all') % default
  cfg.reject = conflist(1:end); % to be removed
else
  cfg.reject = intersect(conflist, cfg.reject); % to be removed
end
fprintf('removing confound %s \n', num2str(cfg.reject));
kprs = setdiff(conflist, cfg.reject); % to be kept
fprintf('keeping confound %s \n', num2str(kprs));

% confound normalization for orthogonality
if ~isfield(cfg, 'normalize') || stcrmp(cfg.normalize, 'yes')
  fprintf('normalizing the confounds, except the constant \n');
  for c = 1:nconf
    SD = std(regr(:,c),0,1);
    if SD == 0
      fprintf('confound %s is a constant \n', num2str(c));
    else
      regr(:,c) = (regr(:,c) - mean(regr(:,c))) / SD;
    end
    clear SD;
  end
elseif stcrmp(cfg.normalize, 'no')
  fprintf('skipping normalization procedure \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM MODEL
%   Y = X * B + err, where Y is data, X is the model, and B are beta's
% which means
%   Best = X\Y ('matrix division', which is similar to B = inv(X)*Y)
% or when presented differently
%   Yest = X * Best
%   Yest = X * X\Y
%   Yclean = Y - Yest (the true 'clean' data is the recorded data 'Y' -
%   the data containing confounds 'Yest')
%   Yclean = Y - X * X\Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isfreq     = ft_datatype(datain, 'freq');
istimelock = ft_datatype(datain, 'timelock');

if istimelock
  switch datain.dimord
    case {'rpt_chan_time', 'subj_chan_time'}
      
      % descriptives
      nrpt  = size(datain.trial, 1);
      nchan = size(datain.trial, 2);
      ntime = size(datain.trial, 3);
      
      % initialize output variable
      dataout       = datain;
      
      if nrpt~=size(regr,1)
        error('the size of your confound matrix does not match with the number of trials/subjects');
      end
      
      % get the data on which the contribution of the confounds has to be estimated
      dat = reshape(datain.trial, [nrpt, nchan*ntime]);
      
      % estimate and remove the confounds
      fprintf('estimating the regression weights and removing the confounds \n');
      beta = regr\dat;                                                        % B = X\Y
      model = regr(:, cfg.reject) * beta(cfg.reject, :);                      % model = confounds * weights = X * X\Y
      Yc = dat - model;                                                       % Yclean = Y - X * X\Y
      
      % put the clean data back into place
      dataout.trial = reshape(Yc, [nrpt, nchan, ntime]); clear Yc;
      
      % update descriptives when already present
      if isfield(dataout, 'var') % remove (old) var
        dataout = rmfield(dataout, 'var');
      end
      if isfield(dataout, 'dof') % remove (old) degrees of freedom
        dataout = rmfield(dataout, 'dof');
      end
      if isfield(dataout, 'avg') % remove (old) avg and reaverage
        fprintf('updating descriptives \n');
        dataout = rmfield(dataout, 'avg');
        tempcfg            = [];
        tempcfg.keeptrials = 'yes';
        dataout = ft_timelockanalysis(tempcfg, dataout); % reaveraging
      end
      
      % make a nested timelock structure that contains the model
      if isfield(cfg, 'model') && strcmp(cfg.model, 'yes')
        fprintf('outputting the model which contains the confounds x weights \n');
        dataout.model.trial   = reshape(model, [nrpt, nchan, ntime]); clear model;
        dataout.model.dimord  = dataout.dimord;
        dataout.model.time    = dataout.time;
        dataout.model.label   = dataout.label;
        if isfield(dataout, 'avg')
          % also average the model
          tempcfg            = [];
          tempcfg.keeptrials = 'yes';
          dataout.model      = ft_timelockanalysis(tempcfg, dataout.model); % reaveraging
        end
      end
      
      % beta statistics
      if isfield(cfg, 'statistics') && strcmp(cfg.statistics, 'yes')
        fprintf('performing statistics on the regression weights \n');
        dfe        = nrpt - nconf;                                              % degrees of freedom
        err        = dat - regr * beta;                                         % err = Y - X * B
        mse        = sum((err).^2)/dfe;                                         % mean squared error
        covar      = diag(regr'*regr)';                                         % regressor covariance
        bvar       = repmat(mse',1,size(covar,2))./repmat(covar,size(mse,2),1); % beta variance
        tval       = (beta'./sqrt(bvar))';                                      % betas -> t-values
        prob       = (1-tcdf(tval,dfe))*2;                                      % p-values
        clear err; clear mse; clear dat; clear regr; clear bvar;
        dataout.stat     = reshape(tval, [nconf, nchan, ntime]); clear tval;
        dataout.prob     = reshape(prob, [nconf, nchan, ntime]); clear prob;
        % FIXME: drop in replace tcdf from the statfun/private dir
      end
      
      % add the beta weights to the output
      dataout.beta     = reshape(beta, [nconf, nchan, ntime]); clear beta;
      
    otherwise
      error('unsupported dimord "%s"', datain.dimord);
  end % switch
  
elseif isfreq
  
  switch datain.dimord
    case {'rpt_chan_freq_time', 'subj_chan_freq_time', 'rpttap_chan_freq_time', 'rpt_chan_freq', 'subj_chan_freq', 'rpttap_chan_freq'}
      
      % descriptives
      nrpt  = size(datain.powspctrm, 1);
      nchan = size(datain.powspctrm, 2);
      nfreq = size(datain.powspctrm, 3);
      ntime = size(datain.powspctrm, 3); % this will be a singleton dimension in case there is no time
      
      % initialize output variable
      dataout       = datain;
      
      if nrpt~=size(confound,1)
        error('the size of your confound matrix does not match with the number of trials/subjects');
      end
      
      % get the data on which the contribution of the confounds has to be estimated
      dat = reshape(datain.powspctrm, [nrpt, nchan*nfreq*ntime]);
      
      % estimate and remove the confounds
      fprintf('estimating the regression weights and removing the confounds \n');
      beta = regr\dat;                                                        % B = X\Y
      Yc   = dat - regr(:, cfg.reject) * beta(cfg.reject, :);                 % Yclean = Y - X * X\Y
      
      % put the clean data back into place
      dataout.powspctrm = reshape(Yc, [nrpt, nchan, nfreq, ntime]); clear Yc;
      dataout.beta      = reshape(beta, [nconf, nchan, nfreq, ntime]);
            
      % beta statistics
      if isfield(cfg, 'statistics') && strcmp(cfg.statistics, 'yes')
      fprintf('performing statistics on the regression weights \n');
      dfe        = nrpt - nconf;                                              % degrees of freedom
      err        = dat - regr * beta;                                         % err = Y - X * B
      mse        = sum((err).^2)/dfe;                                         % mean squared error
      covar      = diag(regr'*regr)';                                         % regressor covariance
      bvar       = repmat(mse',1,size(covar,2))./repmat(covar,size(mse,2),1); % beta variance
      tval       = (beta'./sqrt(bvar))';                                      % betas -> t-values
      prob       = (1-tcdf(tval,dfe))*2;                                      % p-values
      clear err; clear mse; clear bvar; clear dat; clear regr; clear beta;
      dataout.stat  = reshape(tval, [nconf, nchan, nfreq, ntime]); clear tval;
      dataout.prob  = reshape(prob, [nconf, nchan, nfreq, ntime]); clear prob;
      end
      
    otherwise
      error('unsupported dimord "%s"', datain.dimord);
  end % switch
  
  
else
  error('the input data should be either timelock or freq with trials')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath'); % this is helpful for debugging
cfg.version.id   = '$Id: ft_rejectconfound.m 4375 2011-10-07 09:16:45Z arjsto $'; % this will be auto-updated by the revision control system

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();

% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername(); % this is helpful for debugging
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

if hasdata && isfield(datain, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = datain.cfg;
end
clear datain;

% remember the exact configuration details in the output
dataout.cfg = cfg;

% discard the gradiometer information because the weightings have been
% changed
if isfield(dataout, 'grad')
  warning('discarding gradiometer information because the weightings have been changed');
  dataout = rmfield(dataout, 'grad');
end

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', dataout); % use the variable name "data" in the output file
end
