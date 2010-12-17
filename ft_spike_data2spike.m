function [spike] = ft_spike_data2spike(cfg,data)

% FT_SPIKE_DATA2SPIKE converts a continuous spikestation DATA structure to a point
% representation of the data in a spikestation SPIKE structure. 
%
% Configurations options (CFG):
%
%  Configuration options related to selection of spike channel and trials:
%
%   cfg.spikechannel       = string or index of single spike channel to trigger on
%                            (default = 'all'). See FT_CHANNELSELECTION for details

% Martin Vinck 2010 (C). 

if nargin~=2, error('MATLAB:ft_spike_data2spike:data2spike:nargin','Two input arguments required'), end

% put the defaults and check whether there are ununsed fields
defaults.spikechannel   = {'all'};
cfg = ft_spike_sub_defaultcfg(cfg,defaults);

% check if the input data is valid for this function, should be done with CHECKDATA
hasAllFields = all(isfield(data,{'time','trial', 'label'})); 
if ~hasAllFields, error('MATLAB:ft_spike_data2spike:data2spike:wrongInput',...
	'DATA should contain .time, .trial, .label'); 
end

% spike channel selection, avoid check on max number of spikes here, since we may use spike format
% for events as well
cfg.channel = ft_channelselection(cfg.spikechannel, data.label);
spikesel    = match_str(data.label, cfg.channel);
nUnits      = length(spikesel);	% number of units
if nUnits==0, error('MATLAB:ft_spike_data2spike:cfg:spikechannel:noChanSelected',...
      'No spikechannel selected by means of cfg.spikechannel'); 
end

% do the conversion
nTrials 	= length(data.trial);
trialTimes  = zeros(nTrials,2);
for iUnit = 1:nUnits    
    unitIndx = spikesel(iUnit);    
    spikeTimes  = []; % we dont know how large it will be, so use concatenation inside loop
    trialInds   = [];    
    for iTrial = 1:nTrials
        
        % read in the spike times
        [spikeTimesTrial]    = getspiketimes(data, iTrial, unitIndx);        
        nSpikes              = length(spikeTimesTrial);
        spikeTimes           = [spikeTimes; spikeTimesTrial(:)];
        trialInds            = [trialInds; ones(nSpikes,1)*iTrial];
        
        % get the begs and ends of trials
        if iUnit==1, trialTimes(iTrial,:) = data.time{iTrial}([1 end]); end          
    end        
    
    spike.label{iUnit}     = data.label{unitIndx};
    spike.waveform{iUnit}  = [];
    spike.time{iUnit}      = spikeTimes;        
    spike.trial{iUnit}     = trialInds;

    if iUnit==1, spike.trialtime             = trialTimes; end
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
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
spike.cfg     = cfg;


%%%%%%%%%% SUB FUNCTION %%%%%%%%%%
function [spikeTimes spikeIndx] = getspiketimes(data,trial,unit)

% GETSPIKETIMES extracts the spike times and spike samples from a continuous fieldtrip 
% DATA structure in a selected trial for a selected unit.
%
% Inputs:
%   DATA is a contiuous fieldtrip data structure
%   
%   TRIAL is a natural number indicating which trial is selected
%   UNIT  is a natural number indicating which channel is selected
%   
% Outputs:
%   SPIKETIMES contains the spike times, sampled at frequency data.fsample;
%   SPIKEINDX  contains the samples in data.trial{trial}(unit,:) at which we find the
%   spikes.

spikeIndx       = logical(data.trial{trial}(unit,:));            
spikeCount      = data.trial{trial}(unit,spikeIndx);            
spikeTimes      = data.time{trial}(spikeIndx);
multiSpikes     = find(spikeCount>1);                       

% preallocate the additional times that we get from the double spikes
nMultiSpikes = sum(spikeCount(multiSpikes));
[addTimes,addSamples] = deal(zeros(nMultiSpikes,1));

binWidth            = 1/data.fsample; % the width of each bin
halfBinWidth        = binWidth/2; 

% get the additional samples and spike times, we need only loop through the bins
n = 1;
for iBin = multiSpikes(:)' % looping over row vector
    nSpikesInBin = spikeCount(iBin);
    addTimes(n : n+nSpikesInBin-1)   = ones(1,nSpikesInBin)*spikeTimes(iBin);
    addSamples(n : n+nSpikesInBin-1) = ones(1,nSpikesInBin)*spikeIndx(iBin);    
    n = n + nSpikesInBin;
end

% before adding these times, first remove the old ones 
spikeTimes(multiSpikes) = [];
spikeIndx(multiSpikes)  = [];
spikeTimes              = sort([spikeTimes(:); addTimes(:)]);   
spikeIndx               = sort([spikeIndx(:) ; addSamples(:)]);   





        
