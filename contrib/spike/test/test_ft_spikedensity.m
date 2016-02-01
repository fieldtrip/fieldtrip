function test_ft_spikedensity()

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_spikedensity
% TEST ft_spikedensity

data = [];
nTrials = 100;
nSamples = 1001;
R = 0.51;
% create artificial data.
for iTrial = 1:nTrials
    data.trial{iTrial}(1,:) = round(R*rand(1,nSamples));
    data.trial{iTrial}(2,:) = (R*rand(1,nSamples));
end

data.time(1:nTrials) = {linspace(0,1,nSamples)};
data.fsample = 1000;
data.label{1} = 'spk1';
data.label{2} = 'eeg1';
data.hdr = [];
data.cfg.trl = [];

%%
cfgSdf.timwin        = [-0.02 0.02];
cfgSdf.winfunc        = 'gausswin';
cfgSdf.latency       = [0 3];
cfgSdf.keeptrials = 'yes';
cfgSdf.spikechannel = 1
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
%% show that we can also enter a spike input
spike = ft_checkdata(data, 'datatype', 'spike', 'feedback', 'yes');
%%
[sdf2 sdfdata] = ft_spikedensity(cfgSdf,spike);

any(sdf2.trial(:)-sdf.trial(:))>0
%%
nTrials = length(data.trial);
for iTrial = 1 : 5
    figure, plot(sdfdata.time{iTrial}, sdfdata.trial{iTrial}), hold on, 

    spks = find(data.trial{iTrial}(1,:));
    xx = [data.time{iTrial}(1,spks); data.time{iTrial}(1,spks)];
    yy = [ones(1,length(spks))*max(sdfdata.trial{iTrial});ones(1,length(spks))*min(sdfdata.trial{iTrial})];
    hold on
    plot(xx,yy,'r'), 
end
nanmean(sdf.avg) % expect 20 hz
close all
%%
% do the same for rectwin and alphawin
cfgSdf.winfunc = 'rectwin';

tic
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
toc
%%
nTrials = length(data.trial);
for iTrial = 1 : 5
    figure, plot(sdfdata.time{iTrial}, sdfdata.trial{iTrial}), hold on, 

    spks = find(data.trial{iTrial}(1,:));
    xx = [data.time{iTrial}(1,spks); data.time{iTrial}(1,spks)];
    yy = [ones(1,length(spks))*max(sdfdata.trial{iTrial});ones(1,length(spks))*min(sdfdata.trial{iTrial})];
    hold on
    plot(xx,yy,'r'), 
end
close all

% note that the firing rate matches to 1000/41ms ~ 25 hz

%% for alphawin
cfgSdf.winfunc = 'alphawin';

[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
nTrials = length(data.trial);
for iTrial = 1 : 5
    figure, plot(sdfdata.time{iTrial}, sdfdata.trial{iTrial}), hold on, 

    spks = find(data.trial{iTrial}(1,:));
    xx = [data.time{iTrial}(1,spks); data.time{iTrial}(1,spks)];
    yy = [ones(1,length(spks))*max(sdfdata.trial{iTrial});ones(1,length(spks))*min(sdfdata.trial{iTrial})];
    hold on
    plot(xx,yy,'r'), 
end
close all

%% now check if we use different timwinsif we see that back with a gaussian
cfgSdf = [];
cfgSdf.timwin = [0 0.02];
cfgSdf.winfunc = 'gausswin';
cfgSdf.latency = [0 1];

[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
nTrials = length(data.trial);
for iTrial = 1 : 5
    figure, plot(sdfdata.time{iTrial}, sdfdata.trial{iTrial}), hold on, 

    spks = find(data.trial{iTrial}(1,:));
    xx = [data.time{iTrial}(1,spks); data.time{iTrial}(1,spks)];
    yy = [ones(1,length(spks))*max(sdfdata.trial{iTrial});ones(1,length(spks))*min(sdfdata.trial{iTrial})];
    hold on
    plot(xx,yy,'r'), 
end
close all


%% try if this function also works with different latencies as promised in the help

clear
nTrials = 100;
nSamples = 1001;
R = 0.51;
% create artificial data.
data.fsample = 1000;

for iTrial = 1:nTrials
    % first put the latency as an integer of the sample frequency
    latency = 0 + round(20*(rand-0.5))/1000; % is going to include some trials, and exlude some                                      
    latencyEnd = 1 + round(20*(rand-0.5))/1000;
    timeaxis = latency:(1/data.fsample):latencyEnd;
    data.time{iTrial}  = timeaxis;
    data.trial(iTrial) = {round(R*rand(1,length(timeaxis)))};
    latencies(iTrial,:) = minmax(timeaxis);
end
close all

data.label{1} = 'spk1';
data.hdr = [];
data.cfg.trl = [];
%% set the defaults
cfgSdf = [];
cfgSdf.latency = 'minperiod'
%cfgSdf.checkspikechan = 'no'
cfgSdf.timwin = [-0.25 0.25]
cfgSdf.winfunc = 'rectwin'
[sdf sdfdata] = ft_spikedensity(cfgSdf,data); % warning here!
figure, plot(sdf.time,sdf.avg)

%% allow variable trial length, check the dof
cfgSdf = [];
cfgSdf.keeptrials = 'yes';
cfgSdf.timwin = [-0.05 0.05];
cfgSdf.winfunc = 'rectwin';
cfgSdf.latency = 'maxperiod';
cfgSdf.vartriallen = 'yes';
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
figure, plot(sdf.time,sdf.avg)
any(cellfun(@length,data.trial)-cellfun(@length,sdfdata.trial))
% note the discrepancy of one sample, why?
%%
nTrials = length(data.trial);
for iTrial = 1 : 5
    figure, plot(sdfdata.time{iTrial}, sdfdata.trial{iTrial}), hold on, 

    spks = find(data.trial{iTrial});
    xx = [data.time{iTrial}(spks); data.time{iTrial}(spks)];
    yy = [ones(1,length(spks))*max(sdfdata.trial{iTrial});ones(1,length(spks))*min(sdfdata.trial{iTrial})];
    hold on
    plot(xx,yy,'r'), 
end
close all

%sdfdata.trial has maximum length if it fits the window, otherwise it has less length
figure, plot(sdf.time,sdf.avg)
pause(1)
close all

%% do not allow variable trial length
cfgSdf = [];
cfgSdf.timwin = [-0.001 0.001];
cfgSdf.winfunc = 'rectwin';
%cfgSdf.winfuncopt = 0.0;
cfgSdf.latency = [0 1];
cfgSdf.vartriallen = 'yes';
sum(latencies(:,1)>0| latencies(:,2)<1)
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
figure, plot(sdf.time,sdf.avg)
%% do not allow variable trial length
cfgSdf = [];
cfgSdf.timwin = [-0.02 0.02];
cfgSdf.winfunc = 'rectwin';
cfgSdf.latency = [0 1];
cfgSdf.vartriallen = 'no';
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
latencies(sdf.cfg.trials,:)
figure, plot(sdf.time,sdf.avg)
pause(1)
close all


%% check all options one for one
cfgSdf = [];
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
%%
cfgSdf.keeptrials = 'no';
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
cfgSdf = [];
%%
cfgSdf.latency = 'prestim'
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
cfgSdf = [];
%%
cfgSdf.winfunc = 'gauss'
[sdf sdfdata] = ft_spikedensity(cfgSdf,data);
cfgSdf = [];













