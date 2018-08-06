function timelock = fttimelock(this, chanind, timeind, trialind, freqind)
% Method for converting meeg object to Fieldtrip timelock/freq struct
% FORMAT timelock = fttimelock(this, chanind, timeind, trialind, freqind)
%
% The method support both time and TF data and outputs different variants
% of timelock or freq FT struct depending on the dataset type and requested
% data dimensions.
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fttimelock.m 5063 2012-11-16 12:10:41Z vladimir $

if nargin < 2 || isempty(chanind)
    chanind = 1:nchannels(this);
end

if nargin < 3 || isempty(timeind)
    timeind = 1:nsamples(this);
end

if nargin < 4 || isempty(trialind)
    trialind = 1:ntrials(this);
end

if strncmpi(transformtype(this),'TF',2) && ...
        (nargin < 5 || isempty(freqind))
     freqind = 1:nfrequencies(this);
end

timelock             = [];
timelock.label       = chanlabels(this, chanind)';

if isequal(transformtype(this), 'time')
   if isequal(type(this), 'continuous')
            error('For continuous data use ftraw method');
   end
   
   if isequal(type(this), 'single') || length(trialind)>1
       timelock.dimord  = 'rpt_chan_time';
       timelock.trial   =  permute(this.data.y(chanind, timeind, trialind), [3 1 2]);
   else
       timelock.dimord  = 'chan_time';
       timelock.avg     =  spm_squeeze(this.data.y(chanind, timeind, trialind), 3);
   end
   
   timelock.time       = time(this, timeind);
   
elseif strncmpi(transformtype(this),'TF',2)
    if length(timeind)>1
        if isequal(type(this), 'single') || length(trialind)>1
            timelock.dimord    = 'rpt_chan_freq_time';
            timelock.powspctrm = permute(this.data.y(chanind, freqind, timeind, trialind), [4 1 2 3]);
        else
            timelock.dimord     = 'chan_freq_time';
            timelock.powspctrm  =  spm_squeeze(this.data.y(chanind, freqind, timeind, trialind), 3);
        end
        
        timelock.time       = time(this, timeind);
    else
        if isequal(type(this), 'single') || length(trialind)>1
            timelock.dimord    = 'rpt_chan_freq';
            timelock.powspctrm = spm_squeeze(permute(this.data.y(chanind, freqind, timeind, trialind), [4 1 2 3]), 4);
        else
            timelock.dimord     = 'chan_freq';
            timelock.powspctrm  =  spm_squeeze(this.data.y(chanind, freqind, timeind, trialind), [3 4]);
        end
    end
    
    timelock.freq      = frequencies(this, freqind);
else
    error('Unknown transform type.');
end   

if length(trialind)>1
    
    clist      =  condlist(this);
    condlabels = conditions(this, trialind);
    timelock.trialinfo = 0*trialind;
    
    for k = 1:numel(clist)
        fprintf('mapping condition label "%s" to condition code %d\n', clist{k}, k);
        sel = strcmp(clist{k}, condlabels);
        timelock.trialinfo(sel) = k;
    end
    
end

if ~isempty(sensors(this, 'MEG'))
    timelock.grad = sensors(this, 'MEG');
end

if ~isempty(sensors(this, 'EEG'))
    timelock.elec = sensors(this, 'EEG');
end
