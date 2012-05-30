function [spike] = ft_spike_maketrials(cfg,spike)

% FT_SPIKE_MAKETRIALS converts raw timestamps in a SPIKE structure to spike
% times that are relative to an event trigger in an SPIKE structure. This
% is a necessary preprocessing step to use functions such as FT_SPIKE_PSTH.
% The other route is by FT_SPIKE_DATA2SPIKE.
%
% The main function of FT_SPIKE_MAKETRIALS is to create the field
% spike.time and spike.trial, which contain the trial numbers in which the
% spikes were recorded, and the onset and offset of the trial relative to
% trigger time t = 0.
%
% Use as
%   [spike] = ft_spike_maketrials(cfg,spike)
%
% Inputs:
%   The raw spike datatype, obtained from any of the fieldtrip READ functions.
%   Structure spike includes the following spike-format specific fields:
%       - spike.timestamp (1-by-nUnits cell array), raw spike times as from recording,
%         these are not relative to trigger but relative to recording system.
%       - spike.label, cell array containing the labels of spike channels.
%
% Configurations:
%   cfg.trl                     = is an nTrials-by-3 matrix.
%                                Every row contains start (col 1), end (col 2) and offset of the event
%                                trigger in the trial. For example, an offset of -1000 means that the trigger
%                                (t = 0 sec) occurred 1000 timestamps after the trial start.
%   cfg.timestampspersecond     = number of timestaps per second. cfg.timestampspersecond should always
%                                 be explicitly specified.
%
% Outputs:
%   spike.time                  = 1-by-nUnits cell array, containing the spike times in
%                                 seconds relative to the event trigger.
%   spike.trial                 = 1-by-nUnits cella array, containing the trial number for
%                                 every spike telling in which trial it was recorded.
%   spike.trialtime             = nTrials-by-2 matrix specifying the start and end of
%                                 every trial in seconds.
%   spike.trl                   = contains the original matrix of cfg.trl
% Further, reproduced in the output are the original fields spike.timestamp
% and spike.label.

% Copyright (C) 2010, Maritn Vinck; F.C. Donders Centre Nijmegen; University of Amsterdam
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% check if input data is indeed spikeraw format
spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'timestampspersecond', 'trl'});

% ensure that the options are valid
cfg = ft_checkopt(cfg,'timestampspersecond',{'doublescalar'});
cfg = ft_checkopt(cfg,'trl', {'numericvector', 'numericmatrix'});

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
  [ts,indx] = sort(double(ts)); % just sort for safety
  
  % take care of the waveform and make [1 x Samples x Spikes] per default
  ignoreWave = 1;
  if isfield(spike, 'waveform')
    % find the dimension where we have to select
    sz        = size(spike.waveform{iUnit});
    N         = length(ts);
    if length(sz)==2
      if length(N)==sz(2)
        permord = [3 1 2];
      elseif length(N)==sz(1)
        permord = [3 2 1];
      end
      spike.waveform{iUnit} = permute(spike.waveform{iUnit}, permord);
    end
    if ~any(sz==N), 
      ignoreWave = 1;
    else
      ignoreWave = 0;
    end
  end
         
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
    dt = dt/cfg.timestampspersecond + cfg.trl(trialNum,3)/cfg.timestampspersecond;
  else
    dt = [];
  end
  trialDur = double(cfg.trl(:,2)-cfg.trl(:,1))/cfg.timestampspersecond;
  time = [cfg.trl(:,3)/cfg.timestampspersecond (cfg.trl(:,3)/cfg.timestampspersecond + trialDur)]; % make the time-axis
  
  % gather the results
  spike.time{iUnit}   = dt(:)';
  spike.trial{iUnit}  = trialNum(:)';
  spike.trialtime     = time;
  if isfield(spike, 'waveform') && ~ignoreWave
    spike.waveform{iUnit} = spike.waveform{iUnit}(:,:,sel);
    spike.waveformdimord  = 'lead_samples_spikes';
  end
  ts = spike.timestamp{iUnit}(indx(sel));
  spike.timestamp{iUnit} = ts(:)';
end
spike.sampleinfo         = cfg.trl(:,1:2);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous spike
ft_postamble history spike


