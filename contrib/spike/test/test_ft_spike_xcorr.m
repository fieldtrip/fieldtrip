function test_ft_spike_xcorr()

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_ft_spike_xcorr
% TEST ft_spike_xcorr

data = [];
spikesPerTrial = 10;
nTrials = 50;
shapePar = 2;
scalePar = 3;
time       = linspace(0,1,1001);
data.trial(1:nTrials) = {zeros(1,length(time))};
data.time(1:nTrials) = {time};
nUnits = 10;

Spike.timestamp = cell(1,nUnits-1);
Spike.trial = cell(1,nUnits-1);
Spike.time = [];
for iTrial = 1:nTrials    
    iUnit=1
    data.trial{iTrial}(iUnit,:) = rand(1,length(time));

    iUnit = 2;
    spikeTimes = 0.015*gamrnd(shapePar,scalePar,[round(spikesPerTrial*rand) 1]);
    smp = [];
    spikeTimes(spikeTimes>1) = [];
    for iSpike = 1:length(spikeTimes)
      smp(iSpike)        = nearest(time,spikeTimes(iSpike));
    end
    data.trial{iTrial}(iUnit,smp) = 1;       

    iUnit = 3
    spikeTimes = spikeTimes + 0.005;
    smp = [];
    spikeTimes(spikeTimes>1) = [];
    for iSpike = 1:length(spikeTimes)
      smp(iSpike)        = nearest(time,spikeTimes(iSpike));
    end
    data.trial{iTrial}(iUnit,smp) = 1;       
        
    iUnit = 4
    spikeTimes = 0.015*gamrnd(shapePar,scalePar,[round(spikesPerTrial*rand) 1]);
    
    smp = [];
    spikeTimes(spikeTimes>1) = [];
    for iSpike = 1:length(spikeTimes)
      smp(iSpike)        = nearest(time,spikeTimes(iSpike));
    end
    data.trial{iTrial}(iUnit,smp) = 1;       
    
end
data.fsample = 1000;
data.hdr = [];
data.cfg.trl = [];
data.label = {};
for iUnit = 1:4
  data.label{end+1} = strcat('spk',num2str(iUnit));
end
cfgSpike.spikechannel = [2 3 4]
spike = ft_checkdata(data,'datatype', 'spike', 'feedback', 'yes');
%%
% we should find a broad peak now, but also shift predictor with broad peak
cfg = [];
cfg.maxlag   = 0.2;
cfg.keeptrials='yes';
X = ft_spike_xcorr(cfg,spike);
sum(X.xcorr) % adds up to one, correct
figure, plot(X.time,squeeze(X.xcorr(1,2,:))) % should give sharp peak at neg

figure, plot(X.time,squeeze(X.xcorr(1,3,:))) % should give broad peak

figure, plot(X.time,squeeze(X.xcorr(1,1,:))) % should give big peak at center

figure, plot(X.time,squeeze(X.trial(1,1,1,:))) % should give big peak at center
%%
cfg.outputunit = 'raw'
cfg.latency = 'maxperiod'
X = ft_spike_xcorr(cfg,spike);
sum(X.xcorr) % adds up to one, correct


figure, plot(X.time,squeeze(X.xcorr(1,2,:))) % should give sharp peak at neg

figure, plot(X.time,squeeze(X.xcorr(1,3,:))) % should give broad peak

figure, plot(X.time,squeeze(X.xcorr(1,1,:))) % should give big peak at center

figure, plot(X.time,squeeze(X.trial(1,1,2,:))) % should give big peak at center

%%
%cfg.autocorr = 'no';
cfg = [];
cfg.debias = 'yes';
cfg.maxlag   = 0.2;
cfg.method = 'shiftpredictor'
cfg.outputunit = 'raw';
cfg.binsize = 0.001;
cfg.keeptrials = 'yes';
%tic,X = spike_xcorr(cfg,Spike);toc % CrossCorr function gives incorrect results
tic,X2s = ft_spike_xcorr(cfg,spike);toc
cfg.method = 'xcorr'
X2 = ft_spike_xcorr(cfg,spike);toc

%%
figure, plot(squeeze(X2.xcorr(1,3,:)))
hold on
plot(squeeze(X2s.shiftpredictor(1,3,:)),'r')

figure, plot(squeeze(X2.xcorr(1,2,:)))
hold on
plot(squeeze(X2s.shiftpredictor(1,2,:)),'r')

figure, plot(squeeze(X2.xcorr(3,2,:)))
hold on
plot(squeeze(X2s.shiftpredictor(3,2,:)),'r')

