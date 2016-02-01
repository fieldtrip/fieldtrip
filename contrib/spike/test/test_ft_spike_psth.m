function test_ft_spike_psth()

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_spike_psth
% TEST ft_spike_psth

% fixed number of spikes, equal length trials, one LFP channel, gamma distribution psth
% one channel with spikes at known positions
data = [];
spikesPerTrial = 10;
nTrials = 100;
shapePar = 2;
scalePar = 3;
time     = linspace(0,1,1000);
data.trial(1:nTrials) = {zeros(1,length(time))};
data.time(1:nTrials) = {time};
for iTrial = 1:nTrials    
  
  iUnit = 1;
  spikeTimes = [];
  while length(spikeTimes)<spikesPerTrial
    spikeTimes = 0.015*gamrnd(shapePar,scalePar,[spikesPerTrial 1]);
    smp = [];
    spikeTimes(spikeTimes>1) = [];
  end
  spikeTimes = spikeTimes(1:spikesPerTrial);

  for iSpike = 1:length(spikeTimes)
    smp(iSpike)        = nearest(time,spikeTimes(iSpike));
  end  
  data.trial{iTrial}(iUnit,smp) = 1;

  iUnit = 2; % unit with fixed positions
  
  smp = [];
  spikeTimes = linspace(0.1,0.9,9);
  for iSpike = 1:length(spikeTimes)
    smp(iSpike)        = nearest(time,spikeTimes(iSpike));
  end  
  data.trial{iTrial}(iUnit,smp) = 1;
end
data.fsample = 1000;
data.hdr = [];
data.cfg.trl = [];
data.label{1} = 'spk2';
data.label{end+1} = 'spk3';
cfg.spikechannel = 2:3;
spike = ft_checkdata(data,'datatype', 'spike', 'feedback', 'yes');
%%
cfgC.fsample = 1000;
data2 = ft_checkdata(spike,'datatype', 'raw', 'feedback', 'yes', 'fsample', 1000);
data2 = ft_checkdata(data2,'datatype', 'raw', 'feedback', 'yes');
%%
% we compute the psth by calling the psth function
cfgPsth = [];
cfgPsth.binsize       = 0.005;
cfgPsth.outputunit    = 'rate';
psthData = ft_spike_psth(cfgPsth,spike);
%% show that we can also use the data as input
cfgPsth = [];
cfgPsth.binsize       = 0.005;
cfgPsth.outputunit    = 'rate';
psthDataRaw = ft_spike_psth(cfgPsth,data);
psthDataRaw.avg - psthData.avg % should give zeros back

%% make a plot of the psth and compare it to the expected curve
for iChan = 1:2
  figure
  cfgPlot = [];
  cfgPlot.spikechannel = iChan;
  ft_spike_plot_psth(cfgPlot,psthData)
end
figure, plot(time,gampdf(time/0.015,shapePar,scalePar)); % this gives the expected proportion of spikes
pause(1)
close
%%
cfgPlot = [];
cfgPsth.outputunit    = 'spikecount';
psthData = ft_spike_psth(cfgPsth,spike);
figure
cfgPlot.spikechannel = 1;
cfgPlot.ylim = [0 1];
cfgPlot.errorbars = 'conf95%';
ft_spike_plot_psth(cfgPlot,psthData); 
% average is 1, this is correct
%% demonstrate some additional features 
cfgPsth = [];
cfgPsth.binsize = 0.001;
cfgPsth.outputunit  = 'spikecount'
psthData = ft_spike_psth(cfgPsth,spike);

% nothing changes for the second spike chan
figure
cfgPlot.spikechannel = 2;
ft_spike_plot_psth(cfgPlot,psthData); 
% but note how data becomes noisy for the first
figure
cfgPlot.spikechannel = 1;
ft_spike_plot_psth(cfgPlot,psthData); 
%% check with a range of binsizes
for iBinsize = [0 0.001 0.005 0.1 0.5 1 2]
  cfgPsth.binsize = iBinsize;
  try
  psthData = ft_spike_psth(cfgPsth,spike);
  figure
  cfgPlot.spikechannel = 1;
  ft_spike_plot_psth(cfgPlot,psthData); 
  catch
    lasterr
    disp('this should give an error')
  end
end
pause(1)
close
%%
for iBinsize = 1:2
  if iBinsize ==1
    cfgPsth.binsize = 'scott';
  else
    cfgPsth.binsize = 'sqrt';
  end    
  try
  psthData = ft_spike_psth(cfgPsth,spike);
  figure
  cfgPlot.spikechannel = 1;
  ft_spike_plot_psth(cfgPlot,psthData); 
  catch
    lasterr
    disp('this should give an error')
  end
end
pause
close

%% our script should work with variable trial lengths
% create data with variable length at the start

clear
spikesPerTrial = 10;
nTrials   = 100;
shapePar  = 2;
scalePar  = 3;

data = [];
data.time = cell(1, nTrials);
data.trial = cell(1, nTrials);

for iTrial = 1:nTrials    
    
  % create the latency, start and end have a jitter of 100 ms
  latencyBeg = 0 + (100*(rand-0.5))/1000; % is going to include some trials, and exlude some                                      
  latencyEnd = 1 + (100*(rand-0.5))/1000;
  timeAxis = latencyBeg:0.001:latencyEnd;
  n = length(timeAxis);
  
  data.trial{iTrial} = zeros(2, n);

  iUnit = 1;
  spikeTimes = [];
  smp = [];
  while length(spikeTimes)<spikesPerTrial && length(smp)<spikesPerTrial
    spikeTimes = 0.015*gamrnd(shapePar,scalePar,[spikesPerTrial 1]);
  end
  spikeTimes = spikeTimes(1:spikesPerTrial);
  smp = [];
  spikeTimes(spikeTimes>1) = [];
  for iSpike = 1:length(spikeTimes)
    smp(iSpike)        = nearest(timeAxis,spikeTimes(iSpike));    
  end    
  smp(smp>n) = [];
  data.trial{iTrial}(iUnit,smp) = 1;
  
  iUnit = 2;
  spikeTimes = linspace(0,1,10);
  spikeTimes(spikeTimes>timeAxis(end) | spikeTimes<timeAxis(1)) = [];
  smp = [];
  for iSpike = 1:length(spikeTimes)
    smp(iSpike)        = nearest(timeAxis,spikeTimes(iSpike));
  end
  data.trial{iTrial}(iUnit,smp) = 1;
  
  data.time{iTrial}   = timeAxis;
  latencies(iTrial,1) = timeAxis(1);
  latencies(iTrial,2) = timeAxis(end);
end
data.fsample = 1000;
data.hdr = [];
data.cfg.trl = [];
data.label = {};
data.label{end+1} = 'spk2';
data.label{end+1} = 'spk3';
cfg.spikechannel = 1:2;
spike = ft_checkdata(data,'datatype', 'spike', 'feedback', 'yes');

%% test with maximum latency and no variable trial length, dof should be 1 then
% this should give dof of 0 or 1, because latency is maxperiod, so max 1 trial can be
% selected
cfgPsth = [];
cfgPsth.vartriallen  = 'no';
cfgPsth.latency = 'maxperiod'
psthSpike = ft_spike_psth(cfgPsth,spike);
%% accept variable trial length and check
cfgPsth = [];
cfgPsth.vartriallen  = 'yes';
cfgPsth.outputunit = 'spikecount'
cfgPsth.binsize = 0.01;
psthSpike = ft_spike_psth(cfgPsth,spike);
%%
psthSpike.dof % note that dof is decreasing on the borders

%%
figure
cfgPlot = [];
cfgPlot.spikechannel = 2;
h = ft_spike_plot_psth(cfgPlot,psthSpike)
pause(1)
close

%%
% note that the average spike count did not change: this is because if a trial does not
% count in, the dof is also 1 lower
figure
cfgPlot = [];
cfgPlot.spikechannel = 1;
h = ft_spike_plot_psth(cfgPlot,psthSpike)

% note that the side bin has no variance (since dof = 1) and the next bin has high
% sem because dof is lower there
% show how we get nan-padding in psthSpike.trial, and how in trials with right latency we
% have spikes
A = [latencies squeeze(psthSpike.trial(:,2,:))];
A = [[zeros(1,2) psthSpike.time]; A]
pause(1)
close

%% now combine the trial selection with the latency selection
cfgPsth = [];
cfgPsth.trials = find(latencies(:,1)<psthSpike.time(20));
cfgPsth.vartriallen  = 'yes';
cfgPsth.outputunit = 'spikecount';
cfgPsth.latency = [0.2 1]
psthSpike2 = ft_spike_psth(cfgPsth,spike);
% note how the dof is equal now

%%
