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
%     Note that values in cfg.trl get inaccurate above 2^53 (in that case 
%     it is better to use the original uint64 representation)
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init            % this will show the function help if nargin==0 and return an error
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
cfg = ft_checkconfig(cfg, 'allowed', {'datafile', 'dataformat', 'dataset', 'event', 'headerfile', 'headerformat', 'trialfun', 'trlunit', 'timestampspersecond', 'hdr', 'trl'});

if strcmp(cfg.trlunit,'timestamps')
  
  % make a loop through the spike units and make the necessary conversions
  nTrials = size(cfg.trl,1);
  for iUnit = 1:nUnits
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
    
    % take care of the waveform information as well
    hasWave =  isfield(spike, 'waveform') && ~isempty(spike.waveform) && ~isempty(spike.waveform{iUnit});
      
    % make the timestamps relative to the trigger
    trialNum = [];
    sel       = [];
    for iTrial = 1:nTrials
      isVld = find(ts>=trlEvent(iTrial,1) & ts<=trlEvent(iTrial,2));
      if ~isempty(isVld)
        trialNum = [trialNum; iTrial*ones(length(isVld),1)];  %#ok<*AGROW>
      end
      sel   = [sel; isVld(:)]; 
    end

    % subtract the event (t=0) from the timestamps directly
    if ~isempty(trialNum)
      ts = ts(sel);
      dt = double(ts - trlEvent(trialNum,1)); % convert to double only here
      dt = dt/cfg.timestampspersecond + trlDouble(trialNum,3)/cfg.timestampspersecond;
    else
      dt = [];
    end
    trialDur = double(trlEvent(:,2)-trlEvent(:,1))/cfg.timestampspersecond;
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

