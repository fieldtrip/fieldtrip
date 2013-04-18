function [spike] = ft_spike_maketrials(cfg,spike)

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
%
% Inputs:
%     spike: The raw spike datatype, obtained from FT_READ_SPIKE
%
% Configurations:
%
%   cfg.trl  = is an nTrials-by-M matrix, with at least 3 columns:
%     Every row contains start (col 1), end (col 2) and offset of the event
%     trigger in the trial in timestamp or sample units (cfg.trlunit). 
%     For example, an offset of -1000 means that the trigger (t = 0 sec) 
%     occurred 1000 timestamps or samples after the
%     trial start.
%     If more columns are added than 3, these are used to construct the
%     spike.trialinfo field having information about the trial.
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
%     Note that values in cfg.trl get inaccurate above 2^53 (in that case 
%     it is better to use the original uint64 representation)
=======
%     cfg.trl is ideally of the same class as spike.timestamp{} as it avoids round-off
%     errors
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
=======
%     Note:
%     cfg.trl is ideally of the same class as spike.timestamp{} as it avoids round-off
%     errors. In some acquisition systems, the timestamps attain very high values in 
%     uint64 format. If these are represented in a double format, there are round-off
%     errors. As a solution, one should cast the cfg.trl as a uint64 or int64 
%     to avoid the round-off errors.
%     Note that negative numbers are not allowed with uint64. The third column of 
%     cfg.trl may contain a negative offset.
%     To get numerical accuracy one could then cast the cfg.trl as int64.
%     We will then explicitly convert cfg.trl(:,1:2) to uint64 inside the function.
%     
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
=======
>>>>>>> updated comments on ft_spike_maketrials for cfg.trl option
=======
%     Note that values in cfg.trl get inaccurate above 2^53 (in that case 
%     it is better to use the original uint64 representation)
>>>>>>> checking if rounding errors occur in ft_spike_maketrials by seeing whether the cfg.trl exceeds the bitmax
%
%   cfg.trlunit = 'timestamps' (default) or 'samples'. 
%     If 'samples', cfg.trl should 
%     be specified in samples, and cfg.hdr = data.hdr should be specified.
%     This option can be used to reuse a cfg.trl that was used for
%     preprocessing LFP data. 
%     If 'timestamps', cfg.timestampspersecond should be
%     specified, but cfg.hdr should not.
%
%   cfg.hdr     = struct, should be specified if cfg.trlunit = 'samples'.
%     This should be specified as cfg.hdr = data.hdr where data.hdr
%     contains the subfields data.hdr.Fs (sampling frequency of the LFP),
%     data.hdr.FirstTimeStamp, and data.hdr.TimeStampPerSecond.
%
%   cfg.timestampspersecond = number of timestaps per second (for
%     Neuralynx, 1000000 for example). This can be computed for example from
%     the LFP hdr (cfg.timestampspersecond = data.hdr.Fs*data.hdr.TimeStampPerSecond)
%     or is a priori known.
%
% Outputs appended to spike:
%   spike.time                  = 1-by-nUnits cell array, containing the spike times in
%                                 seconds relative to the event trigger.
%   spike.trial                 = 1-by-nUnits cell array, containing the trial number for
%                                 every spike telling in which trial it was recorded.
%   spike.trialtime             = nTrials-by-2 matrix specifying the start and end of
%                                 every trial in seconds.
%   spike.trialinfo             = contains trial information

% Copyright (C) 2010-2013, Martin Vinck
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help            % this will show the function help if nargin==0 and return an error
ft_preamble provenance      % this records the time and memory usage at teh beginning of the function
ft_preamble trackconfig     % this converts the cfg structure in a config object, which tracks the cfg options that are being used
ft_preamble debug           % this allows for displaying or saving the function name and input arguments upon an error

% check if input data is indeed spikeraw format
spike   = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');
nUnits  = length(spike.label);

% process the trlunit field
cfg.trlunit = ft_getopt(cfg,'trlunit','timestamps');
cfg = ft_checkopt(cfg,'trlunit', 'char', {'timestamps', 'samples'});

% process the trl field, which is required
cfg = ft_checkconfig(cfg, 'required', {'trl'});
cfg = ft_checkopt(cfg,'trl', {'numericvector', 'numericmatrix'});
if size(cfg.trl,2)<3,
  warning('cfg.trl should contain at least 3 columns, 1st column start of trial, 2nd column end, 3rd offset, in timestamp or sample units')
end

% check if the cfg.trl is in the right order and whether the trials are overlapping
<<<<<<< HEAD
<<<<<<< HEAD
if strcmp(cfg.trlunit, 'timestamps') && ~all(cfg.trl(:,2)>cfg.trl(:,1))
  warning('the end of some trials does not occur after the beginning of some trials in cfg.trl'); %#ok<*WNTAG>
elseif strcmp(cfg.trlunit, 'samples') && ~all((cfg.trl(:,2))>=cfg.trl(:,1))
  warning('the end of some trials does not occur after the beginning of some trials in cfg.trl'); %#ok<*WNTAG>
end  
  
if size(cfg.trl,1)>1
  if ~all(cfg.trl(2:end,1)>cfg.trl(1:end-1,2))
    warning('your trials are overlapping, trials will not be statistically independent, some spikes will be duplicated'); %#ok<*WNTAG>
  end
end
<<<<<<< HEAD
=======
cfg.trl = cfg.trl;
events  = cfg.trl(:,1:2)'; %2-by-nTrials now
if ~issorted(events(:)), warning('your trials are overlapping, trials will not be statistically independent'); end %#ok<*WNTAG>
if ~issorted(events,'rows'), error('the trials are not in sorted order'); end
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
=======
if ~all(cfg.trl(:,2)>cfg.trl(:,1))
=======
if strcmp(cfg.trlunit, 'timestamps') && ~all(cfg.trl(:,2)>cfg.trl(:,1))
>>>>>>> ensuring correct converstion of samples to timestamps when cfg.trl is in samples for ft_spike_maketrials.m
  warning('the end of some trials does not occur after the beginning of some trials in cfg.trl'); %#ok<*WNTAG>
elseif strcmp(cfg.trlunit, 'samples') && ~all((1+cfg.trl(:,2))>cfg.trl(:,1))
  warning('the end of some trials does not occur after the beginning of some trials in cfg.trl'); %#ok<*WNTAG>
end  
  
if size(cfg.trl,1)>1
  if ~all(cfg.trl(2:end,1)>cfg.trl(1:end-1,2))
    warning('your trials are overlapping, trials will not be statistically independent, some spikes will be duplicated'); %#ok<*WNTAG>
  end
end
>>>>>>> explicitly supporting all possible formats of cfg.trl in ft_spike_maketrials.m

% check if the inputs are congruent: hdr should not be there if unit is timestamps
if strcmp(cfg.trlunit,'timestamps')
  cfg = ft_checkconfig(cfg,'forbidden', 'hdr');
  cfg = ft_checkconfig(cfg, 'required', {'timestampspersecond'});  
  cfg.timestampspersecond = double(cfg.timestampspersecond);  
else
  cfg = ft_checkconfig(cfg, 'required', {'hdr'});      
  if ~isfield(cfg.hdr, 'FirstTimeStamp'), error('cfg.hdr.FirstTimeStamp must be specified'); end
  if ~isfield(cfg.hdr, 'TimeStampPerSample'), error('cfg.hdr.TimeStampPerSample must be specified'); end
  if ~isfield(cfg.hdr, 'Fs'), error('cfg.hdr.Fs, the sampling frequency of the LFP must be specified'); end
end
trlDouble = double(cfg.trl); % this is to compute trial lengths etc.  
%cfg = ft_checkconfig(cfg, 'allowed', {'trlunit', 'timestampspersecond', 'hdr', 'trl'});

if strcmp(cfg.trlunit,'timestamps')
  
  % make a loop through the spike units and make the necessary conversions
  nTrials = size(cfg.trl,1);
  for iUnit = 1:nUnits
<<<<<<< HEAD
<<<<<<< HEAD
    ts = spike.timestamp{iUnit}(:);            
    classTs = class(ts);
    
    % put a warning message if timestamps are doubles but not the right precision
    if (strcmp(classTs, 'double') && any(ts>(2^53))) || (strcmp(classTs, 'single') && any(ts>(2^24)))
      warning('timestamps are of class double but larger than 2^53 or single but larger than 2^24, expecting round-off errors due to precision limitation of doubles');
    end
    
    % check whether trl and ts are of the same class, issue warning if not and it is a problem
    classTrl = class(cfg.trl);
    trlEvent = cfg.trl(:,1:2);
    if ~strcmp(classTs, classTrl)
        flag = 1;
        if strcmp(classTs, 'double') || strcmp(classTrl, 'double')
          mx = 2^53;
          flag = 0;
        end
        if strcmp(classTs, 'single') || strcmp(classTrl, 'single')
          mx = 2^24; % largest precision number
          flag = 0;
        end
        % issue a warning if the class is actually a problem        
        if iUnit==1 && flag==0 && any(cfg.trl(:)>cast(mx, classTrl)) 
          warning('timestamps are of class %s and cfg.trl is of class %s, rounding errors are expected because of high timestamps, converting %s to %s', class(ts), class(cfg.trl), class(cfg.trl), class(ts));
        end
        trlEvent = cast(trlEvent, classTs);
    end
    
=======
    ts = spike.timestamp{iUnit}(:);
<<<<<<< HEAD

>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
=======
=======
    ts = spike.timestamp{iUnit}(:);            
>>>>>>> explicitly supporting all possible formats of cfg.trl in ft_spike_maketrials.m
    classTs = class(ts);
    % put a warning message if timestamps are doubles but not the right precision
    if strcmp(classTs, 'double') && any(ts>(2^53)) || (strcmp(classTs, 'single') && any(ts>(2^24)))
      warning('timestamps are of class double but larger than 2^53 or single but larger than 2^24, expecting round-off errors due to precision limitation of doubles');
    end
    
    classTrl = class(cfg.trl);
    trlEvent = cfg.trl(:,1:2);
    if ~strcmp(classTs, classTrl)
        flag = 1;
        if strcmp(classTs, 'double') || strcmp(classTrl, 'double')
          mx = 2^53;
          flag = 0;
        end
        if strcmp(classTs, 'single') || strcmp(classTrl, 'single')
          mx = 2^24; % largest precision number
          flag = 0;
        end
        if iUnit==1 && flag==0 && any(cfg.trl(:)>cast(mx, classTrl)) 
          % check the maximum to give an indication of the possible error
          warning('timestamps are of class %s and cfg.trl is of class %s, rounding errors are expected because of high timestamps, converting %s to %s', class(ts), class(cfg.trl), class(cfg.trl), class(ts));
        end
        trlEvent = cast(trlEvent, classTs);
    end
    
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
    % take care of the waveform information as well
    hasWave =  isfield(spike, 'waveform') && ~isempty(spike.waveform) && ~isempty(spike.waveform{iUnit});
      
    % make the timestamps relative to the trigger
    trialNum = [];
    sel       = [];
    for iTrial = 1:nTrials
<<<<<<< HEAD
<<<<<<< HEAD
      isVld = find(ts>=trlEvent(iTrial,1) & ts<=trlEvent(iTrial,2));
=======
      if ~strcmp(class(ts), class(cfg.trl))
        isVld = find(double(ts)>=double(events(1,iTrial)) & double(ts)<=double(events(2,iTrial)));
      else
        isVld = find(ts>=events(1,iTrial) & ts<=events(2,iTrial));
      end        
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
=======
      isVld = find(ts>=trlEvent(iTrial,1) & ts<=trlEvent(iTrial,2));
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
      if ~isempty(isVld)
        trialNum = [trialNum; iTrial*ones(length(isVld),1)];  %#ok<*AGROW>
      end
      sel   = [sel; isVld(:)]; 
    end

    % subtract the event (t=0) from the timestamps directly
    if ~isempty(trialNum)
      ts = ts(sel);
<<<<<<< HEAD
<<<<<<< HEAD
      dt = double(ts - trlEvent(trialNum,1)); % convert to double only here
      dt = dt/cfg.timestampspersecond + trlDouble(trialNum,3)/cfg.timestampspersecond;
=======
      if ~strcmp(class(ts), class(cfg.trl))
        warning('timestamps of unit %d are of class %s and cfg.trl is of class %s, rounding errors are possible', iUnit, class(ts), class(cfg.trl));
        dt = double(ts) - double(cfg.trl(trialNum,1));
      else
        dt = double(ts - cfg.trl(trialNum,1)); % convert to double only here
      end
<<<<<<< HEAD
      dt = dt/cfg.timestampspersecond + cfg.trl(trialNum,3)/cfg.timestampspersecond;
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
    else
      dt = [];
    end
    trialDur = double(trlEvent(:,2)-trlEvent(:,1))/cfg.timestampspersecond;
=======
=======
      dt = double(ts - trlEvent(trialNum,1)); % convert to double only here
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
      dt = dt/cfg.timestampspersecond + trlDouble(trialNum,3)/cfg.timestampspersecond;
    else
      dt = [];
    end
<<<<<<< HEAD
    trialDur = double(cfg.trl(:,2)-cfg.trl(:,1))/cfg.timestampspersecond;
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
=======
    trialDur = double(trlEvent(:,2)-trlEvent(:,1))/cfg.timestampspersecond;
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
    time = [trlDouble(:,3)/cfg.timestampspersecond (trlDouble(:,3)/cfg.timestampspersecond + trialDur)]; % make the time-axis

    % gather the results
    spike.time{iUnit}   = dt(:)';
    spike.trial{iUnit}  = trialNum(:)';
    spike.trialtime     = time;
    if hasWave, spike.waveform{iUnit} = spike.waveform{iUnit}(:,:,sel); end
    try spike.unit{iUnit} = spike.unit{iUnit}(sel); end       %#ok<*TRYNC>
    try spike.fourierspctrm{iUnit} = spike.fourierspctrm{iUnit}(sel,:,:); end
    ts = spike.timestamp{iUnit}(sel);
    spike.timestamp{iUnit} = ts(:)';
  end

elseif strcmp(cfg.trlunit,'samples')
   
  nTrials            = size(cfg.trl,1);    
  FirstTimeStamp     = cfg.hdr.FirstTimeStamp;
  TimeStampPerSample = double(cfg.hdr.TimeStampPerSample);
  Fs                 = double(cfg.hdr.Fs);
  cfg.trl            = double(cfg.trl);
  
  [spike.time,spike.trial] = deal(cell(1,nUnits));
  spike.trialtime = zeros(nTrials,2);
  for iUnit = 1:nUnits
    
    % determine the corresponding sample numbers for each timestamp
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
    ts      = spike.timestamp{iUnit}(:);    
    classTs = class(ts);        
    if (strcmp(classTs, 'double') && any(ts>(2^53))) || (strcmp(classTs, 'single') && any(ts>(2^24)))
      warning('timestamps are of class double but larger than 2^53 or single but larger than 2^24, expecting round-off errors due to precision limitation of doubles');
    end
    
    if ~strcmp(classTs, class(FirstTimeStamp))
        flag = 1;
        if strcmp(classTs, 'double') || strcmp(class(FirstTimeStamp), 'double')
          mx = 2^53;
          flag = 0;
        end
        if strcmp(classTs, 'single') || strcmp(class(FirstTimeStamp), 'single')          
          mx = 2^24; % largest precision number
          flag = 0;
        end
        if iUnit==1 && flag==0 && FirstTimeStamp>cast(mx, class(FirstTimeStamp))
           warning('timestamps are of class %s and hdr.FirstTimeStamp is of class %s, rounding errors are possible', class(ts), class(FirstTimeStamp));
        end
        FirstTimeStamp = cast(FirstTimeStamp, classTs);
    end
    sample = double(ts-FirstTimeStamp)/TimeStampPerSample + 1; % no rounding (compare ft_appendspike)
    
    % ensure that cfg.trl is of class double
    if ~strcmp(class(cfg.trl), 'double')
      cfg.trl = double(cfg.trl);
    end
    
    % see which spikes fall into the trials
    waveSel = [];        
=======
    ts = spike.timestamp{iUnit}(:);
    if ~strcmp(class(ts), class(FirstTimeStamp))
      if iUnit==1
        warning('timestamps of unit %d are of class %s and hdr.FirstTimeStamp is of class %s, rounding errors are possible', iUnit, class(ts), class(FirstTimeStamp));
      end
      sample = (double(ts)-double(FirstTimeStamp))/TimeStampPerSample + 1;
    else
      sample = double(ts-FirstTimeStamp)/TimeStampPerSample + 1; % no rounding (compare ft_appendspike)
=======
    ts      = spike.timestamp{iUnit}(:);
=======
    ts      = spike.timestamp{iUnit}(:);    
>>>>>>> explicitly supporting all possible formats of cfg.trl in ft_spike_maketrials.m
    classTs = class(ts);        
    if (strcmp(classTs, 'double') && any(ts>(2^53))) || (strcmp(classTs, 'single') && any(ts>(2^24)))
      warning('timestamps are of class double but larger than 2^53 or single but larger than 2^24, expecting round-off errors due to precision limitation of doubles');
    end
    
    if ~strcmp(classTs, class(FirstTimeStamp))
        flag = 1;
        if strcmp(classTs, 'double') || class(FirstTimeStamp, 'double')
          mx = 2^53;
          flag = 0;
        end
        if strcmp(classTs, 'single') || class(FirstTimeStamp, 'single')          
          mx = 2^24; % largest precision number
          flag = 0;
        end
        if iUnit==1 && flag==0 && FirstTimeStamp>cast(mx, class(FirstTimeStamp))
           warning('timestamps are of class %s and hdr.FirstTimeStamp is of class %s, rounding errors are possible', class(ts), class(FirstTimeStamp));
        end
        FirstTimeStamp = cast(FirstTimeStamp, classTs);
>>>>>>> checking if rounding errors occur in ft_spike_maketrials by seeing whether the cfg.trl exceeds the bitmax
    end
    sample = double(ts-FirstTimeStamp)/TimeStampPerSample + 1; % no rounding (compare ft_appendspike)
<<<<<<< HEAD
    waveSel = [];
>>>>>>> avoiding roundoff errors in ft_spike_maketrials.m
=======
    
    % ensure that cfg.trl is of class double
    if ~strcmp(class(cfg.trl), 'double')
      cfg.trl = double(cfg.trl);
    end
    
    % see which spikes fall into the trials
    waveSel = [];        
>>>>>>> ensuring correct converstion of samples to timestamps when cfg.trl is in samples for ft_spike_maketrials.m
    for iTrial = 1:nTrials
      begsample = cfg.trl(iTrial,1) - 1/2;
      endsample = cfg.trl(iTrial,2) + 1/2;
      sel       = find((sample>=begsample) & (sample<endsample));
      dSample   = sample(sel)-begsample;
      offset    = cfg.trl(iTrial,3)/Fs;               
      tTrial    = dSample/Fs + offset;
      trialNum  = ones(1,length(tTrial))*iTrial;
      trialDur  = (1 + cfg.trl(iTrial,2)-cfg.trl(iTrial,1))/Fs;
      
      spike.time{iUnit}         = [spike.time{iUnit} tTrial(:)'];
      spike.trial{iUnit}        = [spike.trial{iUnit} trialNum];
      if iUnit==1, 
        spike.trialtime(iTrial,:) = [offset offset+trialDur]; 
      end
      waveSel  = [waveSel; sel(:)];
    end     
    
    % select the other fields
    try, spike.waveform{iUnit}      = spike.waveform{iUnit}(:,:,waveSel);      end %#ok<*NOCOM>
    spike.timestamp{iUnit}          = spike.timestamp{iUnit}(waveSel);
    try, spike.unit{iUnit}          = spike.unit{iUnit}(waveSel);              end
    try, spike.fourierspctrm{iUnit} = spike.fourierspctrm{iUnit}(waveSel,:,:); end
  end 
end

if size(cfg.trl,2) > 3
    spike.trialinfo        = double(cfg.trl(:,4:end));
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug            % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig      % this converts the config object back into a struct and can report on the unused fields
ft_postamble provenance       % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and matlab version etc. to the output cfg
ft_postamble previous spike
ft_postamble history spike

