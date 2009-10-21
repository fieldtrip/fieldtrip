function [spike] = spikeanalysis(cfg, data);

% SPIKEANALYSIS performs analysis on spike data
%
% Use as
%   [spike] = spikeanalysis(cfg, data);
%
% The following configuration options are supported
%   cfg.method        = 'rate' (default), or 'spikephase'
%
%   in combination with cfg.method = 'rate',
%     cfg.toi    = the spike-rate is computed in a window surrounding the time-points in cfg.toi
%     cfg.timwin = window width
%     
%     if cfg.toi and cfg.timwin are not specified, the average rate across each trial is computed.
%  
%   in combination with cfg.method = 'spikephase' 
%     cfg.channelcmb    = cell-array, see CHANNELCOMBINATION
%     cfg.bpfilter      = 'no' or 'yes'  bandpass filter
%     cfg.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%     cfg.bpfiltord     = bandpass filter order

% Undocumented local options:
% cfg.bpfilttype
% cfg.channel
% cfg.foi
% cfg.keeptrials
% cfg.output
% cfg.pad
% cfg.previous
% cfg.taper
% cfg.tapsmofrq
% cfg.version

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: spikeanalysis.m,v $
% Revision 1.4  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.3  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.2  2006/01/11 12:47:10  jansch
% implemented option 'rate', which computes either the average rate within each
% trial, or uses freqanalysis_mtmconvol to estimate a windowed spikerate.
%
% Revision 1.1  2005/12/05 10:33:38  roboos
% new implementation after an idea of Thilo
%

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'method'),       cfg.method = 'rate';             end
if ~isfield(cfg, 'channelcmb'),   cfg.channelcmb = {'all', 'all'}; end
if ~isfield(cfg, 'bpfilter'),     cfg.bpfilter = 'yes';            end
if ~isfield(cfg, 'bpfiltord'),    cfg.bpfiltord = 4;               end
if ~isfield(cfg, 'bpfilttype'),   cfg.bpfilttype = 'but';          end
if ~isfield(cfg, 'bpfreq'),       cfg.bpfreq = [30 90];            end

if strcmp(cfg.method, 'rate'),
  ntrials = length(data.trial);
  nchans  = length(data.label);
  spikechan = zeros(nchans,1); 
  for i=1:ntrials
    for j=1:nchans
      spikechan(j) = spikechan(j) + all(data.trial{i}(j,:)==0 | data.trial{i}(j,:)==1);
    end
  end
  chanindx = find(spikechan==ntrials);
  nchans   = length(chanindx);
  label    = data.label(chanindx);

  if ~isfield(cfg, 'toi'),
    %compute the number of spikes within each trial and normalise for the triallength
    rate = zeros(ntrials, nchans);
    for i=1:ntrials
      rate(i, :) = data.fsample*sum(data.trial{i}(chanindx, :), 2)/size(data.trial{i},2)';
    end
  
    spike          = [];
    spike.label    = label(:);
    spike.rate     = rate;
  elseif isfield(cfg, 'toi'),
    if ~isfield(cfg, 'timwin'), error('no timewindow specified'); end;
    %compute the spike-rate based on the DC-bin after fourier-transformation. CAVE: the 
    %data should NOT be baseline corrected.
    tmpcfg           = [];
    tmpcfg.method    = 'mtmconvol';
    tmpcfg.output    = 'pow';
    tmpcfg.keeptrials= 'yes';
    tmpcfg.toi       = cfg.toi;
    tmpcfg.t_ftimwin = cfg.timwin;
    tmpcfg.tapsmofrq = nan;
    tmpcfg.taper     = 'rectwin';
    tmpcfg.foi       = 0;
    tmpcfg.channel   = label;
    tmpcfg.pad       = 'maxperlen';
    tmpfreq          = freqanalysis(tmpcfg, data);
    rate             = sqrt(squeeze(tmpfreq.powspctrm)/2)*data.fsample;    
 
    spike        = [];
    spike.label  = label(:);
    spike.rate   = rate;
  end
elseif strcmp(cfg.method, 'spikephase'),
  % select the combination of spike and lfp channels to be analyzed
  if isfield(cfg, 'channelcmb')
    cfg.channelcmb = channelcombination(cfg.channelcmb, data.label);
  end
  
  % check whether the selection of channels is valid
  ntrials = length(data.trial);
  nchans  = length(data.label);
  spikechan = zeros(nchans,1);
  for i=1:ntrials
    for j=1:nchans
      spikechan(j) = spikechan(j) + all(data.trial{i}(j,:)==0 | data.trial{i}(j,:)==1); 
    end
  end
  spikechan = spikechan==ntrials;
  cmbsel = ones(size(cfg.channelcmb,1), 1);
  for i=1:size(cfg.channelcmb,1)
    cmbsel(i) = spikechan(strmatch(cfg.channelcmb{i,1}, data.label))==1 & ...
                spikechan(strmatch(cfg.channelcmb{i,2}, data.label))==0;
  end
  
  % select only the valid combinations
  cfg.channelcmb = cfg.channelcmb(find(cmbsel),:);
  fprintf('selected %d channel combinations\n', size(cfg.channelcmb,1));
  
  spksel = unique(match_str(data.label, cfg.channelcmb(:,1)));
  lfpsel = unique(match_str(data.label, cfg.channelcmb(:,2)));
  fprintf('selected %d spike channels\n', length(spksel));
  fprintf('selected %d lfp channels\n', length(lfpsel));
  
  % allocate the output structure
  for i=1:size(cfg.channelcmb,1)
    spike.phase{i} = [];
    spike.amplitude{i} = [];
    spike.trlnum{i} = [];
  end
  
  for lfplop=lfpsel(:)'
    fprintf('processing lfp channel %s ', data.label{lfplop});
    for i=1:ntrials
      % perform the bandpass filter and hilbert transform
      lfp = bandpassfilter(data.trial{i}(lfplop,:), data.fsample, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype);
      lfp = hilbert(lfp);
      for spklop=spksel(:)'
        % find the output combination
        outputcmb = find(strcmp(data.label{lfplop}, cfg.channelcmb(:,2)) & strcmp(data.label{spklop}, cfg.channelcmb(:,1)));
        if isempty(outputcmb)
          % this analog channel is not combined with this spike channel
          continue
        else
          fprintf('.');
          % find all spikes in this channel and this trial
          sel = find(data.trial{i}(spklop,:));
          % remember the instantaneous phase and amplitude of the lfp channel
          spike.phase{outputcmb}     = [spike.phase{outputcmb}     phase(lfp(sel))      ];
          spike.amplitude{outputcmb} = [spike.amplitude{outputcmb} abs(lfp(sel))        ]; 
          % remember the trial number in which the spike occurred
          spike.trlnum{outputcmb}    = [spike.trlnum{outputcmb}    i*ones(1,length(sel))];
        end
      end
    end
    fprintf('\n');
  end
  
  % wrap the phase estimate between 0 and 2*pi
  for i=1:size(cfg.channelcmb,1)
    tmp = spike.phase{i};
    tmp = rem(tmp, 2*pi);
    tmp(tmp<0) = tmp(tmp<0) + 2*pi;
    spike.phase{i} = tmp;
  end
  
  % append the other information to the output
  spike.channelcmb = cfg.channelcmb;
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
cfg.version.id = '$Id: spikeanalysis.m,v 1.4 2008/09/22 20:17:44 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
spike.cfg = cfg;

