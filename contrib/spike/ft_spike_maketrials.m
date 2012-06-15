function [spike] = ft_spike_maketrials(cfg,spike, data)

% FT_SPIKE_MAKETRIALS converts raw timestamps in a SPIKE structure to spike
% times that are relative to an event trigger in an SPIKE structure. This
% is a preprocessing step to use functions such as FT_SPIKE_PSTH.
%
% The main function of FT_SPIKE_MAKETRIALS is to create the field
% spike.time and spike.trial, which contain the trial numbers in which the
% spikes were recorded, and the onset and offset of the trial relative to
% trigger time t = 0.
%
% Use as
%   [spike] = ft_spike_maketrials(cfg,spike)
% or 
%   [spike] = ft_spike_maketrials(cfg,spike, data)
%
% Inputs:
%   - spike: The raw spike datatype, obtained from FT_READ_SPIKE
% Optional input:
%   - data: LFP data. 
%     If third input is specified, we create trials based on the trial information present in the LFP. 
%     In this case, data.hdr.FirstTimeStamp and data.hdr.TimeStampPerSample
%     and data.sampleinfo or data.cfg.trl must be present.
%
% Configurations (these do not have to be present if third data input is
% supplied, but are obligatory otherwise.
%
%   cfg.trl  = is an nTrials-by-3 matrix
%     Every row contains start (col 1), end (col 2) and offset of the event
%     trigger in the trial in timestamp units. For example, an offset of -1000 
%     means that the trigger (t = 0 sec) occurred 1000 timestamps after the
%     trial start.
%                               
%   cfg.timestampspersecond = number of timestaps per second (for
%     Neuralynx, 1000000 for example). This can be computed for example from
%     the LFP hdr (cfg.timestampspersecond = data.hdr.Fs*data.hdr.TimeStampPerSecond)
%
% Outputs appended to spike:
%   spike.time                  = 1-by-nUnits cell array, containing the spike times in
%                                 seconds relative to the event trigger.
%   spike.trial                 = 1-by-nUnits cella array, containing the trial number for
%                                 every spike telling in which trial it was recorded.
%   spike.trialtime             = nTrials-by-2 matrix specifying the start and end of
%                                 every trial in seconds.
%   spike.trl                   = contains the original matrix of cfg.trl
%   spike.sampleinfo            = beginning and end of trials in timestamps

% Copyright (C) 2010-2012, Martin Vinck
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
nUnits  = length(spike.label);
  
if nargin==2
  % ensure that the required options are present
  cfg = ft_checkconfig(cfg, 'required', {'timestampspersecond', 'trl'});

  % ensure that the options are valid
  cfg = ft_checkopt(cfg,'timestampspersecond',{'doublescalar'});
  cfg = ft_checkopt(cfg,'trl', {'numericvector', 'numericmatrix'});

  % make sure that the cfg.trl indeed has three columns
  if size(cfg.trl,2)~=3,
    error('cfg.trl should contain 3 columns, 1st column start of trial, 2nd column end, 3rd offset, in timestamp units')
  end

  % make a loop through the spike units and make the necessary conversions
  nTrials = size(cfg.trl,1);
  cfg.trl = double(cfg.trl);
  for iUnit = 1:nUnits
    ts = double(spike.timestamp{iUnit}(:));

    % take care of the waveform and make [1 x Samples x Spikes] per default
    hasWave =  isfield(spike, 'waveform') && ~isempty(spike.waveform) && ~isempty(spike.waveform{iUnit});

    % check if the events are overlapping or not
    events = cfg.trl(:,1:2)'; %2-by-nTrials now
    if ~issorted(events(:))
      warning('Your trials are overlapping, trials will not be statistically independent');
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

    % subtract the event (t=0) from the timestamps directly
    if ~isempty(trialNum)
      ts  	 = ts(sel);
      dt = ts - cfg.trl(trialNum,1); % error if empty
      dt = dt/cfg.timestampspersecond + cfg.trl(trialNum,3)/cfg.timestampspersecond;
    else
      dt = [];
    end
    trialDur = (cfg.trl(:,2)-cfg.trl(:,1))/cfg.timestampspersecond;
    time = [cfg.trl(:,3)/cfg.timestampspersecond (cfg.trl(:,3)/cfg.timestampspersecond + trialDur)]; % make the time-axis

    % gather the results
    spike.time{iUnit}   = dt(:)';
    spike.trial{iUnit}  = trialNum(:)';
    spike.trialtime     = time;
    if hasWave, spike.waveform{iUnit} = spike.waveform{iUnit}(:,:,sel); end
    try spike.unit{iUnit} = spike.unit{iUnit}(sel); end      
    try spike.fourierspctrm{iUnit} = spike.fourierspctrm{iUnit}(sel,:,:); end
    ts = spike.timestamp{iUnit}(indx(sel));
    spike.timestamp{iUnit} = ts(:)';
  end
  spike.sampleinfo         = cfg.trl(:,1:2);

  % do the general cleanup and bookkeeping at the end of the function
  ft_postamble trackconfig
  ft_postamble callinfo
  ft_postamble previous spike
  ft_postamble history spike
else
   
  data = ft_checkdata(data,'type', 'raw', 'feedback', 'yes');
  try
    trl = double(ft_findcfg(data.cfg, 'trl'));
  catch
    try  
      trl = double(data.sampleinfo);
    end
  end
  if isempty(trl), error('could not find the trial information in the continuous data'); end
  nTrials = size(trl,1);
  
  try
    FirstTimeStamp     = double(data.hdr.FirstTimeStamp);
    TimeStampPerSample = double(data.hdr.TimeStampPerSample);
  catch
    error('could not find the timestamp information in the continuous data');
  end
  spike.sampleinfo = double(trl(:,1:2)-1)*TimeStampPerSample + FirstTimeStamp;
  
  [spike.time,spike.trial] = deal(cell(1,nUnits));
  spike.trialtime = zeros(nTrials,2);
  for iUnit = 1:nUnits
    
    % determine the corresponding sample numbers for each timestamp
    ts = double(spike.timestamp{iUnit}(:));
    sample = (ts-FirstTimeStamp)/TimeStampPerSample + 1; % no rounding (compare ft_appendspike)
    waveSel = [];
    for iTrial = 1:nTrials
      begsample = trl(iTrial,1);
      endsample = trl(iTrial,2);
      sel       = find((sample>=begsample) & (sample<=endsample));
      dSample   = sample(sel)-begsample;
      tTrial    = dSample/data.fsample + data.time{iTrial}(1);
      trialNum  = ones(1,length(tTrial))*iTrial;
                   
      spike.time{iUnit}         = [spike.time{iUnit} tTrial(:)'];
      spike.trial{iUnit}        = [spike.trial{iUnit} trialNum];
      if iUnit==1, 
        spike.trialtime(iTrial,:) = [data.time{iTrial}(1) data.time{iTrial}(end)]; 
      end
      finalsel = sel;
      waveSel  = [waveSel; finalsel(:)];
      % construct the sample info based on timestamps
    end 
    % gather the results
    
    try, spike.waveform{iUnit} = spike.waveform{iUnit}(:,:,waveSel); end
    spike.timestamp{iUnit} = spike.timestamp{iUnit}(waveSel);
    try, spike.unit{iUnit}      = spike.unit{iUnit}(waveSel); end
    try, spike.fourierspctrm{iUnit} = spike.fourierspctrm{iUnit}(waveSel,:,:); end
  end 
  ft_postamble trackconfig
  ft_postamble callinfo
  ft_postamble previous spike
  ft_postamble history spike
end


