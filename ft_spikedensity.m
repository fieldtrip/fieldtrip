function [sdf sdfdata] = density(cfg,data)

% FT_SPIKESTATION_DENSITY computes the spikedensityfunction of the spiketrains
% by convolving with a window specified by CFG.WINDOW.
%
% Use as
%    [sdf sdfdata] = ft_spikestation_density(cfg, data)
% or 
%    [sdf]        = ft_spikestation_density(cfg, data)
%    In this case only the average and variance of spikedensityfunction across trials will
%    be computed and individual trials are not kept (SDFDATA has those).
%
%   See cfg.winfunc below for more information on the specific use.
%
% Inputs:
%   DATA should be organised in a structure as obtained from the APPENDSPIKE or 
%   FT_SPIKESTATION_SPIKE2DATA or PREPROCESSING function. 
%
%   Configurations (CFG):
%
%   cfg.timwin         = [begin end], time of the smoothing kernel (default = [-0.1 0.1])
%                        If cfg.winfunc = @alphawin, cfg.timwin(1) will be set to 0.
%   cfg.outputunit     = 'rate' (default) or 'spikecount'. This determines the physical unit
%                        of our spikedensityfunction, either in firing rate or in
%                        spikecount.
%   cfg.winfunc        = (a) string or function handle, type of window to convolve with (def = @gauss).
%                            Options should be set with cfg.winfuncopt
%                        - @gauss (default)
%                        - @alphawin, given by win = x*exp(-x/timeconstant)
%                        - For standard window functions in the signal processing toolbox see
%                          WINDOW.
%                        (b) vector of length nSamples, used directly as window 
%   cfg.winfuncopt     = options that go with cfg.winfunc
%                        For cfg.winopt = @alpha: the timeconstant in seconds (default = 0.005s)
%                        For cfg.winopt = @gauss: the standard devision in seconds (default =
%                                         1/4 of window duration in seconds)
%                        For cfg.winopt = @wname with @wname any standard window function
%                                         see window opts in that function and add as cell array                       
%   cfg.latency        = [begin end] in seconds, 'maxperiod' (default), 'minperiod', 
%                        'prestim'(t>=0), or 'poststim' (t>=0).
%   cfg.vartriallen    = 'yes' (default) or 'no'. 
%                        'yes' - accept variable trial lengths and use all available trials
%                         and the samples in every trial. Missing values will be ignored in
%                         the computation of the average and the variance.
%                        'no'  - only select those trials that fully cover the window as
%                         specified by cfg.latency.
%   cfg.trials         =  numeric or logical selection of trials (default = 'all')
%   cfg.keeptrials     = 'yes' or 'no' (default). If 'yes', we store the trials in a matrix
%                         in output SDF as well.
% Outputs:
%   - SDF is a structure similar to TIMELOCK (output from TIMELOCKANALYSIS) and can be used
%     in TIMELOCKSTATISTICS for example.
%   - SDFDATA is a raw DATA type structure that can be used itself in all functions that support
%     raw data input (such as TIMELOCKANALYSIS, FREQANALYSIS).

% Copyright (C) 2010, Martin Vinck
% TODO: check that SDFDATA is indeed completely compatible!


% set the defaults
defaults.timwin         = {[-0.05 0.05]};          
defaults.trials         = {'all'};           
defaults.latency        = {'maxperiod'};                
defaults.spikechannel   = {'all'};           
defaults.winfunc        = {'gauss'};    
defaults.winfuncopt     = {[]};
defaults.vartriallen    = {'yes' 'no'};            
defaults.outputunit     = {'rate' 'spikecount'};          
defaults.keeptrials     = {'no'  'yes'};
cfg		 = ft_spike_sub_defaultcfg(cfg,defaults);

% check if the input data is valid for this function
hasAllFields = all(isfield(data,{'time','trial', 'label','fsample'})); 
if ~hasAllFields, error('MATLAB:spikestation:density:wrongInput',...
	'DATA should contain .time, .trial, .label and .fsample fields'); 
end

% select the units
cfg.channel = ft_channelselection(cfg.spikechannel, data.label);
spikesel    = match_str(data.label, cfg.channel);
nUnits      = length(spikesel); % number of spike channels
if nUnits==0, error('MATLAB:spikestation:density:cfg:spikechannel:noSpikeChanSelected',...
                    'No spikechannel selected by means of cfg.spikechannel'); 
end

% get the number of trials or change DATA according to cfg.trials
if  strcmp(cfg.trials,'all')    
    cfg.trials = 1:length(data.trial);
elseif islogical(cfg.trials)
    cfg.trials = find(cfg.trials); 
elseif ~isrealvec(cfg.trials);
    error('MATLAB:spikestation:density:cfg:trials:unknownOption',...
    'cfg.trials should be logical or numerical selection or string "all"');  
end
cfg.trials = sort(cfg.trials(:));
if max(cfg.trials)>length(data.trial),error('MATLAB:spikestation:density:cfg:trials:maxExceeded',...
   'maximum trial number in cfg.trials should not exceed length of DATA.trial')
end
if isempty(cfg.trials), 
  errors('MATLAB:spikestation:density:cfg:trials','No trials were selected in cfg.trials'); 
end

% determine the duration of each trial
begTrialLatency = cellfun(@min,data.time(cfg.trials));
endTrialLatency = cellfun(@max,data.time(cfg.trials));

% select the latencies
if strcmp(cfg.latency,'minperiod')
      cfg.latency = [max(begTrialLatency) min(endTrialLatency)];
elseif strcmp(cfg.latency,'maxperiod')
      cfg.latency = [min(begTrialLatency) max(endTrialLatency)];
elseif strcmp(cfg.latency,'prestim')
      cfg.latency = [min(begTrialLatency) 0];
elseif strcmp(cfg.latency,'poststim')
      cfg.latency = [0 max(endTrialLatency)];
elseif ~isrealvec(cfg.latency)||length(cfg.latency)~=2
  error('MATLAB:spikestation:density:cfg:latency',...
    'cfg.latency should be "max", "min", "prestim", "poststim" or 1-by-2 numerical vector');
end
if ~isrealvec(cfg.timwin)||length(cfg.timwin)~=2
  error('MATLAB:spikestation:density:cfg:timwin',...
    'cfg.latency should be 1-by-2 numerical vector');
end
if cfg.latency(1)>=cfg.latency(2),
   error('MATLAB:spikestation:density:cfg:latency:wrongOrder',...
   'cfg.latency should be a vector in ascending order, i.e., cfg.latency(2)>cfg.latency(1)');
end
if (cfg.latency(1) < min(begTrialLatency)), cfg.latency(1) = min(begTrialLatency); 
  warning('MATLAB:spikestation:density:begLatencyTooEarly',...
          'Correcting begin latency of averaging window');
end
if (cfg.latency(2) > max(endTrialLatency)), cfg.latency(2) = max(endTrialLatency);
  warning('MATLAB:spikestation:density:endLatencyTooLate',...
          'Correcting end latency of averaging window');
end

% start processing the window information
if strcmp(cfg.winfunc,'alphawin') % now force start of window to be positive.
    warning('MATLAB:spikestation:density:cfg:timwin:alphawin:timwinNeg',...
      'cfg.timwin(1) should be a non-negative number if cfg.winfunc = @alphawin')
    cfg.timwin(1) = 0;  
end

% construct the time-axis for the window
if cfg.timwin(1)>0 || cfg.timwin(2)<0 || cfg.timwin(1)>=cfg.timwin(2)
  error('MATLAB:spikestation:density:cfg:timwin:noInclusionTimeZero',...
        'Please specify cfg.timwin(1)<=0 and cfg.timwin(2)>=0 and cfg.timwin(2)>cfg.timwin(1)')
end
sampleTime    = 1/data.fsample;
winTime       = [fliplr(0:-sampleTime:cfg.timwin(1)) (sampleTime):sampleTime:cfg.timwin(2)];
nLeftSamples  = length(find(winTime<0));  
nRightSamples = length(find(winTime>0));  
nSamplesWin   = length(winTime);
if nSamplesWin==1, warning('MATLAB:spikestation:density:cfg:timwin:winLengthOne',...
  'Number of samples in selected window is exactly one, so no smoothing applied')
end

% construct the window
% construct the window
if strcmp(cfg.winfunc,'gauss') 
   if  isempty(cfg.winfuncopt), cfg.winfuncopt{1} = 0.25*diff(cfg.timwin); end
   win = exp(-(winTime.^2)/(2*cfg.winfuncopt{1}^2));
elseif strcmp(cfg.winfunc,'alphawin')
   if isempty(cfg.winfuncopt),  cfg.winfuncopt{1} = 0.005; end    
   win = winTime.*exp(-winTime/cfg.winfuncopt{1});   
elseif ischar(cfg.winfunc) || strcmp(class(cfg.winfunc),'function_handle') 
   if isempty(cfg.winfuncopt)
     win = feval(cfg.winfunc,nSamplesWin); 
   else
     win = feval(cfg.winfunc,nSamplesWin, cfg.winfuncopt{:}); 
   end
elseif isnumeric(cfg.winfunc) % only do error check here
  if  ~isrealvec(cfg.winfunc) || length(cfg.winfunc)~=nSamplesWin
       error('MATLAB:spikestation:density:cfg:window:wrongSize', '%s\n%s%d',...
      'cfg.winfunc should be 1-by-N vector, with N equal to the number', ...
      'of samples as determined by cfg.timwin, namely',length(cfg.winfunc))
  end     
else   error('MATLAB:spikestation:density:cfg:window:unknownOption','%s\n%s',...
      'cfg.winfunc should be "gausswin_private", "alphawin", window function (string or handle)',...
      'or 1-by-N vector');
end
win    = win(:).'./sum(win);    % normalize the window to 1
winDur = max(winTime) - min(winTime); % duration of the window

% check which trials will be used based on the latency 
% at this point cfg.trials has selected trials and begTrialLatency is accordingly
overlaps      = endTrialLatency>=(cfg.latency(1)+winDur) & begTrialLatency<=(cfg.latency(2)-winDur);
if strcmp(cfg.vartriallen,'no') % only select trials that fully cover our latency window
    startsLater    = single(begTrialLatency>single(cfg.latency(1)));
    endsEarlier    = single(endTrialLatency<single(cfg.latency(2)));
    hasWindow      = ~(startsLater | endsEarlier); % it should not start later or end earlier
else
    hasWindow      = true(length(cfg.trials),1); % in case vartriallen = "yes"
end
trialSel          = overlaps(:) & hasWindow(:); % trials from cfg.trials we select further
cfg.trials        = cfg.trials(trialSel);       % cut down further on cfg.trials
begTrialLatency   = begTrialLatency(trialSel);  % on this variable as well
nTrials           = length(cfg.trials);         % the actual number of trials we will use
if isempty(cfg.trials),warning('MATLAB:ft_spikedensity:cfg:trials:noneSelected',...
   'no trials were selected, please check cfg.trials')
end

% calculates the samples we are shifted wrt latency(1)
samplesShift      = zeros(1,nTrials); 
sel               = (begTrialLatency-cfg.latency(1))>0; % otherwise 0 samples to be padded
samplesShift(sel) = round(data.fsample*(begTrialLatency(sel)-cfg.latency(1))); 

% create the time axis for the spike density
time           = cfg.latency(1):(1/data.fsample):cfg.latency(2); 
cfg.latency(2) = time(end);   % this is the used latency that should be stored in cfg again
maxNumSamples  = length(time);

% preallocate the sum, squared sum and degrees of freedom
[s,ss]   = deal(zeros(nUnits, maxNumSamples)); % sum and sum of squares
dof      = zeros(1, length(s)); 

if (strcmp(cfg.keeptrials,'yes')), singleTrials = zeros(nTrials,nUnits,size(s,2)); end

% preallocate depending on whether nargout is 1 or 2
if nargout==2, [sdfdata.trial(1:nTrials) sdfdata.time(1:nTrials)] = deal({[]}); end
for iTrial = 1:nTrials    
    origTrial  = cfg.trials(iTrial);   % this is the original trial we use for DATA input
    timeAxis   = data.time{origTrial}; % get the time axis for this trial
    sampleSel  = nearest(timeAxis, cfg.latency(1)) : nearest(timeAxis, cfg.latency(2)); 
    nSamples   = length(sampleSel);
    trialTime  = timeAxis(sampleSel);  % select the relevant portion of time

    for iUnit = 1:nUnits        
        unitIndx = spikesel(iUnit); % index in data.label                
        dat      = data.trial{origTrial}(unitIndx,sampleSel); % get the data
        if any(dat)
          y = conv(full(dat),win);  % use convolution to get the raw spike density
        else
          y = zeros(1,nSamples + nSamplesWin - 1); % no spikes; no convolution needed
        end
        
        if strcmp(cfg.outputunit, 'rate') 
          y = y*data.fsample;  % normalize to the sampling rate, to get it in Hz.
        else
          y = y*nSamplesWin;   % now maximum becomes N (length window)
        end                
        y = y(nLeftSamples+1 : end-nRightSamples); % delete additional points we get with conv                
        y([1:nLeftSamples end-nRightSamples+1:end]) = NaN; % assign NaN at borders
        
        if nargout==2           
           sdfdata.trial{iTrial}(iUnit,:) = y;
           sdfdata.time{iTrial}(iUnit,:)  = trialTime; % write back original time axis
        end
        
        dofsel = true(1,length(y));
        if strcmp(cfg.vartriallen,'yes')     
            padLeft  = zeros(1, samplesShift(iTrial));
            padRight = zeros(1,(maxNumSamples - nSamples - samplesShift(iTrial)));
            ySingleTrial = [NaN(size(padLeft)) y NaN(size(padRight))];
            y        = [padLeft y padRight];            
            dofsel   = logical([padLeft dofsel padRight]);
        end
        s(iUnit,:)        = s(iUnit,:)  + y;            % compute the sum
        ss(iUnit,:)       = ss(iUnit,:) + y.^2;         % compute the squared sum
        
        % count the number of samples that went into the sum
        dof(dofsel) = dof(dofsel) + 1;        
        
        if strcmp(cfg.keeptrials,'yes'), singleTrials(iTrial,iUnit,:) = ySingleTrial; end
    end    
    % remove the trial from data in order to avoid buildup in memory
    data.trial{origTrial} = [];
    data.time{origTrial}  = [];
end

dofMat            = dof(ones(nUnits,1),:);

% give back a similar structure as timelockanalysis
sdf.avg                  = s ./ dofMat; 
sdf.var                  = (ss - s.^2./dofMat)./(dofMat-1); % sumPsth.^2 ./ dof = dof .* (sumPsth/dof).^2
sdf.dof                  = dofMat;
sdf.time                 = time;
sdf.fsample         	 = data.fsample; 
sdf.label(1:nUnits)      = data.label(spikesel);
if (strcmp(cfg.keeptrials,'yes'))
  sdf.trial = singleTrials;
  sdf.dimord = 'rpt_chan_time';
else
  sdf.dimord = 'chan_time';
end

% create a new structure that is a standard raw data spikestation structure itself
if nargout>1
    sdfdata.fsample              = data.fsample;
    sdfdata.label(1:nUnits)      = data.label(spikesel);
    sdfdata.hdr                  = data.hdr;
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
% remember the configuration details of the input data
if isfield(data,'cfg'), cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
sdf.cfg     = cfg;
if nargout>1
  sdfdata.cfg = cfg;
end

