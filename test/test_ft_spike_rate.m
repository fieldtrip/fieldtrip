function test_ft_spike_rate()

% TEST test_ft_spike_rate

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
spike = ft_spike_data2spike(cfg,data);

%%
cfgRate = [];
cfgRate.outputunit = 'spikecount';
RateCnt = ft_spike_rate(cfgRate,spike);


RateCnt.avg
RateCnt.trial