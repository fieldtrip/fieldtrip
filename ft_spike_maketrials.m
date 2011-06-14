function [spike] = ft_spike_maketrials(cfg,spike)

% FT_SPIKE_MAKETRIALS converts raw timestamps in a SPIKE structure to spike times
% that are relative to an event trigger in an SPIKE structure. 
% This is a necessary preprocessing step to use functions such as FT_SPIKE_PSTH.
% The other route is by FT_SPIKE_DATA2SPIKE.
%
% Inputs:
%   The raw spike datatype, obtained from any of the fieldtrip READ functions.
%   Structure spike includes the following spike-format specific fields:
%       - spike.timestamp (1-by-nUnits cell array), raw spike times as from recording,
%         these are not relative to trigger but relative to recording system.
%       - spike.label, cell array containing the labels of spike channels.
%
%   The main function of FT_SPIKE_MAKETRIALS is to create the field spike.time and
%   spike.trial, which contain the trial numbers in which the spikes were recorded, and
%   the onset and offset of the trial relative to trigger time t = 0.   
%
% Configurations:
%   cfg.trl                     = is an nTrials-by-3 matrix.
%                                Every row contains start (col 1), end (col 2) and offset of the event 
%                                trigger in the trial. For example, an offset of -1 sec means that trigger 
%                                (t = 0 sec) occurred 1 second after the trial start.
%   cfg.timestampspersecond     = number of timestaps per second. cfg.timestampspersecond should always 
%                                 be explicitly specified.
%
%   Outputs:
%   spike.time                  = 1-by-nUnits cell array, containing the spike times in
%                                 seconds relative to the event trigger.
%   spike.trial                 = 1-by-nUnits cella array, containing the trial number for
%                                 every spike telling in which trial it was recorded.
%   spike.trialtime             = nTrials-by-2 matrix specifying the start and end of
%                                 every trial in seconds.
%   spike.trl                   = contains the original matrix of cfg.trl
%   Further, reproduced in spike are the original fields spike.timestamp and spike.label

%   Martin Vinck (C) 2010, F.C. Donders Centre Nijmegen, University of Amsterdam

if nargin<2, 
  error('MATLAB:ft_spike_maketrials:nargin','I can not work with less than 2 inputs'), 
end

% detect whether the format of spike is correct
hasAllFields = all(isfield(spike, {'timestamp' 'label'}));
if ~hasAllFields, error('MATLAB:ft_spike_maketrials:wrongStructInput',...
      'HELP. Input SPIKE should be struct with .timestamp and .label fields.')
end

% check whether all are of right format
correctInp = iscell(spike.timestamp) & iscell(spike.label);
if ~correctInp, error('MATLAB:ft_spike_maketrials:wrongStructInput',...
    'HELP. .timestamp and .label should be cell arrays.')
end

% make sure that the user explicitly specifies the timestamps per second
if ~isfield(cfg,'timestampspersecond') 
    error('MATLAB:ft_spike_maketrials:cfg:timestampspersecond',...
            'How many timestamps are there in a second?')  
end

% make sure that the user explicitly specifies the timestamps per second
if ~isfield(cfg,'trl') 
    error('MATLAB:ft_spike_maketrials:cfg:timestampspersecond',...
            'Please give in a cfg.trl structure, I have no clue when your monkey dances')  
end


% make sure that the cfg.trl indeed has three columns
if size(cfg.trl,2)~=3, 
   error('MATLAB:ft_spike_maketrials:TRL',...
         'HELP. TRL should contain 3 columns, 1st column start of trial, 2nd column end, 3rd offset')
end

% make a loop through the spike units and make the necessary conversions
nTrials = size(cfg.trl,1);
nUnits  = length(spike.label);
cfg.trl = double(cfg.trl);
for iUnit = 1:nUnits
      
    ts = spike.timestamp{iUnit}(:);
    ts = sort(double(ts)); % just sort for safety
    
    % check if the events are overlapping or not
    events = double(cfg.trl(:,1:2))'; %2-by-nTrials now
    if ~issorted(events(:))
        warning('MATLAB:ft_spike_maketrials:trialoverlap',...
        'Your trials are overlapping, trials will not be statistically independent'); 
    end        

    % make the timestamps relative, use different algorithm when overlapping (fast & slow)
    trialNum = [];
    sel       = [];
    for iTrial = 1:nTrials
      isVld = find(ts>=events(1,iTrial) &ts<=events(2,iTrial));
      if ~isempty(isVld)
          trialNum = [trialNum; iTrial*ones(length(isVld),1)];
      end
      sel   = [sel; isVld(:)];
    end
   
    % subtract the event (t=0) from the timestamps directly, this can be double or uint64
    if ~isempty(trialNum)  
        ts  	 = ts(sel);        
        dt = ts - cfg.trl(trialNum,1); % error if empty
        dt = dt/cfg.timestampspersecond + cfg.trl(trialNum,3);    
    else
        dt = [];
    end
    trialDur = double(cfg.trl(:,2)-cfg.trl(:,1))/cfg.timestampspersecond;
    time = [cfg.trl(:,3) (cfg.trl(:,3) + trialDur)]; % make the time-axis                
    
    % gather the results
    spike.time{iUnit}   = dt;
    spike.trial{iUnit}  = trialNum;
    spike.trialtime     = time;  
end
spike.trl = cfg.trl;

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
cfg.previous = [];
try, cfg.previous = spike.cfg; end
% remember the exact configuration details in the output 
spike.cfg = cfg;


