function [wave,spike] = ft_spike_waveform_descriptive(cfg,spike)

% FT_SPIKE_WAVEFORM_DESCRIPTIVE computes descriptive parameters on
% waveform (mean and variance), and performs operations like reallignment, outlier rejection,
% invertation, normalization and interpolation (see configurations).
%
% Use as
%   [wave] = ft_spike_waveform_descriptive(cfg, spike)
% Or
%   [wave, spike] = ft_spike_waveform_descriptive(cfg, spike)
% The input SPIKE should be organised as the SPIKE datatype (see FT_DATATYPE_SPIKE)
%
% Configurations:
%   cfg.rejectoutliers   = 'yes' (default) or 'no': takes away waveforms with too late peak, and no
%                           rising AP towards peak of other waveforms
%   cfg.normalize        = 'yes' (default) or 'no': normalizes all waveforms
%                           to have peak-to-through amp of 2
%   cfg.interp           = 'yes' (default) or 'no'. If 'yes', we interpolate
%   cfg.interpfreq       = frequency of interp resampling. For example, 10
%                          (default)
%                          creates 10 subsamples per sample.
%   cfg.allign           = 'yes' (def). or 'no'. If 'yes', we allign all waves to
%                          maximum
%   cfg.fsample          = sampling frequency of waveform time-axis.
%                          Obligatory field.
%   cfg.spikechannel     = See FT_CHANNELSELECTION for details.
%   cfg.keeptrials       = 'yes' or 'no' (default). If 'yes', we keep the
%                          changed waveform in output
%
% Outputs:
%  Wave.avg   = average waveform
%  Wave.time  = time of waveform axis
%  Wave.var   = variance of waveform
%  Wave.dof   = number of spikes contributing to average
%
% Spike structure if two outputs are desired: waveform is replaced by interpolated and
% cleaned waveforms, removing also their associated time-stamps and data.
%
%  Copyright (C) 2012, Martin Vinck & Thilo Womelsdorf
%

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'fsample'});

% get the default options
cfg.allign          = ft_getopt(cfg, 'allign','yes');
cfg.interp          = ft_getopt(cfg, 'interp','yes');
cfg.interpfreq      = ft_getopt(cfg, 'interpfreq',10);
cfg.spikechannel    = ft_getopt(cfg, 'spikechannel', 'all');
cfg.rejectonpeak    = ft_getopt(cfg,'rejectonpeak', 'yes');
cfg.rejectclippedspikes  = ft_getopt(cfg,'rejectclippedspikes', 'yes');
cfg.normalize       = ft_getopt(cfg,'normalize', 'yes');
cfg.keeptrials      = ft_getopt(cfg,'keeptrials', 'yes');

% ensure that the options are valid
cfg = ft_checkopt(cfg,'allign','char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'rejectclippedspikes','char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'rejectonpeak','char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'interp','char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'normalize','char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg,'keeptrials', 'char', {'yes', 'no'});
cfg = ft_checkopt(cfg,'fsample', 'double');

spike = ft_checkdata(spike,'datatype', 'spike', 'feedback', 'yes');

cfg.channel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel    = match_str(spike.label, cfg.channel);
nUnits      = length(spikesel);
if nUnits==0, error('No spikechannel selected by means of cfg.spikechannel'); 
end

[mnWaveform, varWaveform, dofWaveform] = deal([]);
for iUnit = 1:nUnits
  fprintf('Processing waveforms for the %d th unit \n', iUnit);
  % check if we should get the first dimension first
  spikeindx = spikesel(iUnit);
  waves     = spike.waveform{spikeindx};  
  [nLeads nSamples nSpikes] = size(waves);
  samples = 0:nSamples-1;
  
  % detect if we need to invert the waveform or not
  maxSample = round(2*nSamples/3);
  mn        = nanmean(nanmean(waves,3),1); % nLeads by nSamples
  [vl,iup]   = nanmax(mn(1:maxSample));
  [vl,idown] = nanmin(mn(1:maxSample));
  if iup>idown, waves = -waves; end % if maximum follows after minimum, invert.
  
  maxSample = round(2*nSamples/3);
  mn        = nanmean(nanmean(waves,3),1); % nLeads by nSamples
  [vl,iup]   = nanmax(mn(1:maxSample));
  [vl,idown] = nanmin(mn(1:maxSample));

  % discard waveforms where derivative is not positive until peak index
  if strcmp(cfg.rejectonpeak,'yes') && idown>iup   
    fprintf('Removing units with strange rise and late peak\n')
    % reject the ones that do not have a rising potential to the peak index
    % do this for all four leads at the same time
    mnOverLead = nanmean(waves,1);
    d = find(squeeze(nansum(diff(mnOverLead(:,1:iup,:),[],2),2)));    
    rm1 = find(d<0);
    
    % these have a later max than min: reject, this removes late peaks.
    mnOverLead = nanmean(waves,1);        
    [vl,iu]   = nanmax(mnOverLead,[],2);
    [vl,id]   = nanmin(mnOverLead,[],2);
    %waves(:,:,find(iu>id)) = [];
    rm2 = find(iu>id);
  else
    rm1 = [];
    rm2 = [];
  end
  
  if strcmp(cfg.rejectclippedspikes,'yes')
    fprintf('Removing spikes whose APs were clipped\n')
        
    % reject clipped spikes
    dl = [];
    for iWave = 1:nSpikes
      mnOverLead = nanmean(waves(:,:,iWave),1);            
      df         = diff(mnOverLead);
      idx        = find(df==0)+1;
      if ~isempty(idx), 
        if all(nanmax(mnOverLead(idx))>=mnOverLead) | all(nanmin(mnOverLead(idx))<=mnOverLead)
          dl = [dl iWave];
        end
      end
    end
    rm3 = dl;
  else
    rm3 = [];
  end    
  toRemove = unique([rm1(:); rm2(:); rm3(:)]);
  waves(:,:,toRemove) = [];    
  fprintf('Removing %d spikes from unit %s\n', length(toRemove), spike.label{spikeindx});
  [nLeads nSamples nSpikes] = size(waves);          
  
  % allign the waveforms automatically to the peak index
  % the same allignment must be done for the four leads of a trode
  if strcmp(cfg.allign,'yes')
    if strcmp(cfg.rejectonpeak,'no')
      warning('alligning without rejecting outliers is dangerous');
    end
    mnWave             = nanmean(waves,1);
    [ignore,alignIndx] = nanmax(mnWave(:,:,:),[],2); % maximum across units
    padLen       = nSamples+1 ;
    samplesShift = -padLen:padLen;
    wavesShift   = NaN(nLeads,length(samplesShift),nSpikes);
    peakIndx     = nearest(samplesShift,0);
    startIndx    = peakIndx - iup + 1;
    actIndx      = iup-alignIndx+startIndx;
    for j = 1:nSpikes
        wavesShift(:,actIndx(j):actIndx(j)+nSamples-1,j) = squeeze(waves(:,:,j));
    end
    % --- get x axis right before upscaling and parameter extraction
    time  = samplesShift/cfg.fsample;
    waves = wavesShift;
  else
    time = samples/cfg.fsample;
  end
  
  % make it outlier sensitive by actually looking at the modus  
  if strcmp(cfg.interp,'yes')
     t2 = linspace(time(1), time(end), length(time)*cfg.interpfreq);
     Wn = zeros(nLeads,length(t2),nSpikes);
     warning off
     for k=1:nSpikes
       for iLead = 1:nLeads
         Wn(iLead,:,k) = interp1(time,squeeze(waves(iLead,:,k)),t2);
       end
     end
     warning on
     waves = Wn;
     time = t2;
  end
  
  dof = sum(~isnan(waves),3);
  sd  = nanstd(waves,[],3);  
  %maxSample  = round(2*size(waves,2)/3);
  mn         = nanmean(nanmean(waves,3),1); % nLeads by nSamples
  [vl,iup]   = nanmax(mn(1:end));
  [vl,idown] = nanmin(mn(1:end));
  mn         = nanmean(waves,3); % nLeads by nSamples
  [nLeads nSamples nSpikes] = size(waves);          

  % --- normalize amplitude ratio of spike waveforms if requested
  if strcmp(cfg.normalize,'yes')
      r = mn(:,iup)-mn(:,idown);
      mn = 2*mn./ repmat(r,1,nSamples); %makes it have amp of 2
      mn = mn + -repmat(mn(:,idown)+1,1,nSamples); % put the minus on -1
      sd = sd*2./repmat(r,1,nSamples); % since std(cX) = c std(X);
  end
   
  % --- compute the average waveform here and put in a structure for the next function
  mnWaveform(iUnit,:,:)  = mn;
  varWaveform(iUnit,:,:) = sd.^2;
  dofWaveform(iUnit,:,:) = dof;
  if strcmp(cfg.keeptrials,'yes')
    wave.waveformdimord = '{chan}_lead_spike_time';
    wave.waves{iUnit}   = waves;
  end
  if nargout==2
    spike.waveform{spikeindx} = waves;
    spike.waveformtime = time;
    try, spike.timestamp{spikeindx}(toRemove) = [];end    
    try, 
      if isfield(spike,'trial'), spike.trial{spikeindx}(toRemove) = []; end, 
    end    
    try, 
      if isfield(spike,'unit'), spike.unit{spikeindx}(toRemove) = []; end, 
    end    
        
    try, 
      if isfield(spike,'time'), spike.time{spikeindx}(toRemove) = []; end,     
    end    
    try, 
      if isfield(spike,'fourierspctrm'), spike.fourierspctrm{spikeindx}(toRemove,:,:) = []; end,     
    end    
  end
end
    
wave.time  = time;
wave.avg   = mnWaveform;
wave.dof   = dofWaveform;
wave.var   = varWaveform;
wave.label = spike.label(spikesel);
wave.dimord = 'chan_lead_time';

ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous spike
ft_postamble history spike wave