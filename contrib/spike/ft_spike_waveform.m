function [wave, spike] = ft_spike_waveform(cfg, spike)

% FT_SPIKE_WAVEFORM computes descriptive parameters on
% waveform (mean and variance), and performs operations like realignment, outlier rejection,
% invertation, normalization and interpolation (see configurations).
%
% Use as
%   [wave] = ft_spike_waveform(cfg, spike)
% Or
%   [wave, spike] = ft_spike_waveform(cfg, spike)
% The input SPIKE should be organised as the SPIKE datatype (see FT_DATATYPE_SPIKE)
%
% Configurations:
%   cfg.rejectonpeak     = 'yes' (default) or 'no': takes away waveforms with too late peak, and no
%                           rising AP towards peak of other waveforms
%   cfg.rejectclippedspikes = 'yes' (default) or 'no': removes spikes that
%                           saturated the voltage range. 
%   cfg.normalize        = 'yes' (default) or 'no': normalizes all
%   waveforms
%                           to have peak-to-through amp of 2
%   cfg.interpolate      = double integer (default = 1). Increaes the
%                          density of samples by a factor cfg.interpolate
%   cfg.align            = 'yes' (def). or 'no'. If 'yes', we align all waves to
%                          maximum
%   cfg.fsample          = sampling frequency of waveform time-axis.
%                          Obligatory field.
%   cfg.spikechannel     = See FT_CHANNELSELECTION for details.
%
% Outputs:
%   Wave.avg   = average waveform
%   Wave.time  = time of waveform axis
%   Wave.var   = variance of waveform
%   Wave.dof   = number of spikes contributing to average
%
% Spike structure if two outputs are desired: waveform is replaced by interpolated and
% cleaned waveforms, removing also their associated time-stamps and data.

%  Copyright (C) 2012, Martin Vinck & Thilo Womelsdorf
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
ft_preamble init
ft_preamble provenance spike
ft_preamble trackconfig

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', {'fsample'});

% support the typo in this cfg option that was present in older versions of this function
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1814
cfg = ft_checkconfig(cfg, 'renamed', {'allign', 'align'});

% get the default options
cfg.align               = ft_getopt(cfg, 'align','yes');
cfg.interpolate         = ft_getopt(cfg, 'interpolate', 1);
cfg.spikechannel        = ft_getopt(cfg, 'spikechannel', 'all');
cfg.rejectonpeak        = ft_getopt(cfg,'rejectonpeak', 'yes');
cfg.rejectclippedspikes = ft_getopt(cfg,'rejectclippedspikes', 'yes');
cfg.normalize           = ft_getopt(cfg,'normalize', 'yes');

% ensure that the options are valid
cfg = ft_checkopt(cfg, 'align','char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'rejectclippedspikes','char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'rejectonpeak','char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'interpolate','doublescalar');
cfg = ft_checkopt(cfg, 'normalize','char', {'yes', 'no'});
cfg = ft_checkopt(cfg, 'spikechannel',{'cell', 'char', 'double'});
cfg = ft_checkopt(cfg, 'fsample', 'double');

cfg = ft_checkconfig(cfg, 'allowed', {'align', 'rejectclippedspikes', 'rejectonpeak', 'interpolate', 'normalize', 'spikechannel', 'fsample'});

spike = ft_checkdata(spike, 'datatype', 'spike', 'feedback', 'yes');

cfg.channel = ft_channelselection(cfg.spikechannel, spike.label);
spikesel    = match_str(spike.label, cfg.channel);
nUnits      = length(spikesel);
if nUnits==0, error('No spikechannel selected by means of cfg.spikechannel'); end

[mnWaveform, varWaveform, dofWaveform] = deal([]);
for iUnit = 1:nUnits
  fprintf('Processing waveforms for the unit %d \n', iUnit);
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
    
    fprintf('Removing spikes with strange rise and late peak\n')
    
    % reject the ones that do not have a rising potential to the peak index
    mnOverLead = nanmean(waves,1);% do this for all four leads at the same time    
    d = squeeze(nansum(diff(mnOverLead(:,1:iup,:),[],2),2));    
    rm1 = find(d<0);
    
    % these have a later max than min: reject, this removes late peaks.
    mnOverLead = nanmean(waves,1);        
    [vl,iu]   = nanmax(mnOverLead,[],2);
    [vl,id]   = nanmin(mnOverLead,[],2);
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
        if all(nanmax(mnOverLead(idx))>=mnOverLead) || all(nanmin(mnOverLead(idx))<=mnOverLead)
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
  fprintf('Keeping %d spikes from unit %s\n', nSpikes, spike.label{spikeindx});
  
  % align the waveforms automatically to the peak index
  % the same alignment must be done for the four leads of a trode
  if strcmp(cfg.align,'yes')
    if strcmp(cfg.rejectonpeak,'no')
      warning('aligning without rejecting outliers is dangerous');
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
    
    % --- get time axis right
    time  = samplesShift/cfg.fsample;
    waves = wavesShift;
  else
    time = samples/cfg.fsample;
  end
  
  % interpolate waveforms 
  if cfg.interpolate > 1
     t2 = linspace(time(1), time(end), round(length(time)*cfg.interpolate));
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
  mnWaveform(iUnit,1:nLeads,:)  = mn;
  varWaveform(iUnit,1:nLeads,:) = sd.^2;
  dofWaveform(iUnit,1:nLeads,:) = dof;
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

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble previous   spike
ft_postamble provenance wave spike
ft_postamble history    wave spike

