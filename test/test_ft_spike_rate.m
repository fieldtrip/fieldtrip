function test_ft_spike_rate()

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_spike_rate
% TEST ft_spike_rate

clear
spikesPerTrial = 1;
nTrials = 100;
shapePar = 2;
scalePar = 3;
data = [];
time       = linspace(0,1,1000);
data.trial(1:nTrials) = {zeros(1,length(time))};
data.time(1:nTrials) = {time};
for iTrial = 1:nTrials    
  for iUnit = 1:4
      data.trial{iTrial}(iUnit,:) = rand(1,1000)<iTrial/100;
      if iUnit==4
        data.trial{iTrial}(iUnit,:) = rand(1,1000);
      end
  end
end
data.fsample = 1000;
data.hdr = [];
data.cfg.trl = [];
data.label{1} = 'chan1';
data.label{end+1} = 'chan2';
data.label{end+1} = 'chan3';
data.label{end+1} = 'chan4';

% show that the psth works also with the poisson format
cfg.spikechannel = 1:3;
spike = ft_checkdata(data,'datatype', 'spike', 'feedback', 'yes');

%%
cfgRate = [];
cfgRate.outputunit = 'spikecount';
cfgRate.keeptrials = 'yes';
cfgRate.vartriallen = 'no';
cfgRate.latency = [-10 20];
RateCnt1 = ft_spike_rate(cfgRate,spike);

%%
cfgRate = [];
cfgRate.outputunit = 'spikecount';
cfgRate.keeptrials = 'no';
cfgRate.trials = 1:2:30;
RateCnt1 = ft_spike_rate(cfgRate,spike);



cfgRate.trials = 32:2:60;
RateCnt2 = ft_spike_rate(cfgRate,spike);
cfgRate.trials = 61:2:99;
RateCnt3 = ft_spike_rate(cfgRate,spike);

cfg.stimuli = [0 pi/8 pi/4];
cfg.method = 'orientation';
stat = ft_spike_rate_orituning(cfg,RateCnt1,RateCnt2,RateCnt3);
%%

cfgRate = [];
cfgRate.outputunit = 'spikecount';
cfgRate.keeptrials = 'yes';
cfgRate.trials = 1:5:30;
RateCnt1 = ft_spike_rate(cfgRate,spike);
disp('expect about 10 times the number of spikes as the trial number')
[RateCnt1.cfg.trials(:) RateCnt1.trial]
