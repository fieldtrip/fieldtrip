function test_ft_spike_jpsth()

% MEM 2gb
% WALLTIME 00:10:00

% TEST test_ft_spike_jpsth
% TEST ft_spike_jpsth

%% 
% 1 channelcmb with same latency, gamma
% 1 channelcmb with fixed spike times
% 1 channelcmb with shifted latency to see effect
% 1 lfp channel
% latencies equal
spikesPerTrial = 10;
nTrials        = 10;
shapePar       = 2;
scalePar       = 3;
time           = linspace(0,1,1000);
data.trial(1:nTrials) = {zeros(1,length(time))};
data.time(1:nTrials) = {time};
for iTrial = 1:nTrials    
  iUnit = 1; % lfp channel
  data.trial{iTrial}(iUnit,:) = rand(1,1000);
  
  iUnit = 2;
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
  data.trial{iTrial}(3,smp) = 1;

  iUnit = 4; % unit with fixed positions
  spikeTimes = spikeTimes + 0.01;
  smp = [];
  for iSpike = 1:length(spikeTimes)
    smp(iSpike)        = nearest(time,spikeTimes(iSpike));
  end  
  data.trial{iTrial}(iUnit,smp) = 1;

  iUnit = 5
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

        
  iUnit = 6;
    smp = [];
  spikeTimes = linspace(0.1,0.9,9);
  for iSpike = 1:length(spikeTimes)
    smp(iSpike)        = nearest(time,spikeTimes(iSpike));
  end  
  data.trial{iTrial}(iUnit:iUnit+1,smp) = 1;
  
  iUnit = 8;
    smp = [];
  spikeTimes = linspace(0.1,0.9,9)+0.01;
  for iSpike = 1:length(spikeTimes)
    smp(iSpike)        = nearest(time,spikeTimes(iSpike));
  end  
  data.trial{iTrial}(iUnit:iUnit+1,smp) = 1;
  
  
  
end
data.fsample = 1000;
data.hdr = [];
data.cfg.trl = [];
for iChan = 1:9
data.label{iChan} = strcat('chan',num2str(iChan));
end

data.fsample = 1000;
data.hdr = [];
data.cfg.trl = [];
%%

% show that the psth works also with the poisson format
cfgDE.spikechannel = 2:9;
spike = checkdata(data,'datatype', 'spike', 'feedback', 'yes')


%%
% we compute the psth by calling the psth function
cfgPsth = [];
cfgPsth.binsize       = 0.003; 
cfgPsth.outputunit    = 'rate';
psth = ft_spike_psth(cfgPsth,spike);
%%
cfg = [];
cfg.method = 'jpsth';
tic,jpsth = ft_spike_jpsth(cfg,psth);toc
%%
figure
cfg = [];
cfg.interpolate = 4;
cfg.window = 'gausswin';
cfg.winlen = 0.02
cfg.channelcmb = {jpsth.label{1},jpsth.label{2}};
ft_spike_plot_jpsth(cfg,jpsth)
pause(1)
close all
%%
figure
cfg.channelcmb = {jpsth.label{1},jpsth.label{3}};
ft_spike_plot_jpsth(cfg,jpsth)
% note the shift along the diagonal
figure
cfg = [];
cfg.channelcmb  = {'chan7' 'chan6'};
ft_spike_plot_jpsth(cfg,jpsth)

figure
cfg = [];
cfg.channelcmb  = {'chan8' 'chan6'};
ft_spike_plot_jpsth(cfg,jpsth)
pause(1)
close all

%% test with the sdf
cfg = [];
cfg.keeptrials = 'yes';
cfg.timwin = [-0.01 0.01];
[sdf] = ft_spikedensity(cfg,data);
cfg = [];
cfg.channelcmb = {'chan2' 'chan3'};
tic,jpsth = ft_spike_jpsth(cfg,sdf);toc
% note that it works on SDF as well!
figure
cfg = [];
cfg.channelcmb = {'chan3' 'chan2'};
ft_spike_plot_jpsth(cfg,jpsth)
pause(1)
close all


%%
cfg = [];
cfg.method = 'shiftpredictor';
tic,jpsth = ft_spike_jpsth(cfg,psth);toc
cfg = [];
figure
cfg.channelcmb = {'chan5' 'chan2'};
ft_spike_plot_jpsth(cfg,jpsth)
pause(1)
close all

% note the randomness around 0, shiftpredictor gives same result
%%
figure
cfg.channelcmb = {'chan3' 'chan2'};
ft_spike_plot_jpsth(cfg,jpsth)

pause(1)
close all

 %% test the normalization option
cfg = [];
cfg.normalization = 'yes';
tic,jpsth = ft_spike_jpsth(cfg,psth);toc
cfg = [];
for iCmb1 = 1:7
  for iCmb2 = iCmb1:8
    figure
    cfg.channelcmb = {jpsth.label{iCmb1},jpsth.label{iCmb2}};
    ft_spike_plot_jpsth(cfg,jpsth)
  end
end
figure
cfg.channelcmb = {'chan8', 'chan9'};
ft_spike_plot_jpsth(cfg,jpsth)

nanmax(jpsth.jpsth(:))
nanmin(jpsth.jpsth(:))
pause(1)
close all

% note how the normalized lies between -1 and 1
%%
cfg = [];
cfg.normalization = 'yes';
cfg.method = 'shiftpredictor';
tic,jpsth = ft_spike_jpsth(cfg,psth);toc
cfg = [];
%%
for iCmb1 = 1:7
  for iCmb2 = iCmb1:8
    figure
    cfg.channelcmb = {jpsth.label{iCmb1},jpsth.label{iCmb2}};
    ft_spike_plot_jpsth(cfg,jpsth)
  end
end
pause(1)
close all

%%
figure
cfg.channelcmb = {'chan8' ,'chan9'};
ft_spike_plot_jpsth(cfg,jpsth)

figure
cfg.channelcmb = {'chan3' ,'chan2'};
ft_spike_plot_jpsth(cfg,jpsth)
pause(1)
close all

% note how the normalized lies between -1 and 1

%% check how the jpsth behaves if psth was based on variable trial length
% create data with variable length at the start

clear
spikesPerTrial  = 10;
nTrials         = 100;
shapePar        = 2;
scalePar        = 3;

data = [];
data.time = cell(1, nTrials);
data.trial = cell(1, nTrials);

for iTrial = 1:nTrials    
    
  % create the latency, start and end have a jitter of 100 ms
  latencyBeg = 0 + (100*(rand-0.5))/1000; % is going to include some trials, and exlude some                                      
  latencyEnd = 1 + (100*(rand-0.5))/1000;
  timeAxis   = latencyBeg:0.001:latencyEnd;
  n          = length(timeAxis);
  
  data.trial{iTrial} = zeros(3, n);

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
  spikeTimes = spikeTimes + 0.02;
  smp = [];
  spikeTimes(spikeTimes>1) = [];
  for iSpike = 1:length(spikeTimes)
    smp(iSpike)        = nearest(timeAxis,spikeTimes(iSpike));    
  end    
  smp(smp>n) = [];
  data.trial{iTrial}(iUnit,smp) = 1;
  
  iUnit = 3;
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
data.label{end+1} = 'spk4';
spike = ft_checkdata(data, 'datatype', 'spike', 'feedback', 'yes');

%%
% we compute the psth by calling the psth function
cfgPsth = [];
cfgPsth.binsize       = 0.003; 
cfgPsth.outputunit    = 'rate';
psth = ft_spike_psth(cfgPsth,spike);
%%
psth.dof
squeeze(psth.trial(2,:,:))

%%
cfg = [];
cfg.trials = 1:10;
tic,jpsth = ft_spike_jpsth(cfg,psth);toc

%% check whether plotting still works fine with nans
figure
cfg = [];
cfg.channelcmb = {jpsth.label{1},jpsth.label{2}};
ft_spike_plot_jpsth(cfg,jpsth)
pause(1)
close all

