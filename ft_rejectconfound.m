function dataout = ft_rejectconfound(cfg, datain)

% FT_REGRESSCONFOUND estimates the regression weight of a set of confounds
% and removes the estimated contribution from the single-trial data
%
% Use as
%  timelock = ft_rejectconfound(cfg, timelock)
% or
%  freq = ft_rejectconfound(cfg, freq)
% where datain comes from FT_TIMELOCKANALYSIS or FT_FREQANALYSIS with
% keeptrials = 'yes' and cfg is a configuratioun structure that should
% contain
%   cfg.confound    = matrix, Ntrials X Nconfounds
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

% Copyrights (C) 2011, Robert Oostenveld
%
% $Id$

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

% get the options
confound = ft_getopt(cfg, 'confound');  % there is no default value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isfreq     = ft_datatype(datain, 'freq');
istimelock = ft_datatype(datain, 'timelock');

if istimelock
  switch datain.dimord
    case {'rpt_chan_time', 'subj_chan_time'}
      nrpt  = size(datain.trial, 1);
      nchan = size(datain.trial, 2);
      ntime = size(datain.trial, 3);
      
      if nrpt~=size(confound,1)
        error('the size of your confound matrix does not match with the number of trials/subjects');
      end
      
      % get the data on which the contribution of the confounds has to be estimated
      dat = reshape(datain.trial, [nrpt, nchan*ntime]);
      
      % the model is
      %   Y = X*B + err
      % which means
      %   Best = X\Y
      %   Yest = X * Best
      %   Yclean = Y - Yest = Y - X * X\Y
      
      % estimate and remove the confounds
      fprintf('estimating the regression weight and removing the confounds \n');
      dat = dat - confound * (confound \ dat);
      
      if false
        % FIXME, the definition of beta and B is not clear to me
        B = confound \ dat;
        dataout.B = B;
      end
      
      % put the clean data back into place
      dataout = datain;
      dataout.trial = reshape(dat, [nrpt, nchan, ntime]);
      
    otherwise
      error('unsupported dimord "%s"', datain.dimord);
  end % switch
  
elseif isfreq
  
  switch datain.dimord
    case {'rpt_chan_freq_time', 'subj_chan_freq_time', 'rpttap_chan_freq_time', 'rpt_chan_freq', 'subj_chan_freq', 'rpttap_chan_freq'}
      nrpt  = size(datain.powspctrm, 1);
      nchan = size(datain.powspctrm, 2);
      nfreq = size(datain.powspctrm, 3);
      ntime = size(datain.powspctrm, 4); % this will be a singleton dimension in case there is no time
      
      if nrpt~=size(confound,1)
        error('the size of your confound matrix does not match with the number of trials/subjects');
      end
      
      % get the data on which the contribution of the confounds has to be estimated
      dat = reshape(datain.powspctrm, [nrpt, nchan*nfreq*ntime]);
      
      % estimate and remove the confounds
      fprintf('estimating the regression weight and removing the confounds \n');
      dat = dat - confound * (confound \ dat);
      
      % put the clean data back into place
      dataout = datain;
      dataout.powspctrm = reshape(dat, [nrpt, nchan, nfreq, ntime]);
      
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
cfg.version.id   = '$Id$'; % this will be auto-updated by the revision control system

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

% remember the exact configuration details in the output
dataout.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', dataout); % use the variable name "data" in the output file
end
