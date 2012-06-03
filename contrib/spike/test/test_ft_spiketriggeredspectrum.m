function test_ft_spiketriggeredspectrum()

% TEST test_ft_spiketriggeredspectrum
% ft_spiketriggeredspectrum

%% create the data, in which spikes fall at known phases, and in which our spike times
% fall precicely on the samples. In this case there is no phase error
data = [];
label{1} = 'lfp1';
label{2} = 'lfp2';
label{3} = 'spk';
data.label = label;
cfg.timwin = [-0.1 0.1];
cfg.foilim = [10 90];
data.fsample = 1000;

% time: use multiple trials, each from 0 to 2 seconds
% fs = 1000;
NTRIALS = 20;
time(1:NTRIALS) = {(0:0.001:2)};

% enter a sinusoid with a known phase, without noise for now
SINFREQ = 50; %hz
trial(1:NTRIALS) = {[rand(1,length(time{1}))+cos(SINFREQ*2*pi*time{1});rand(1,length(time{1}))+sin(SINFREQ*2*pi*time{1})]}; % so a spike at t = 0 should mean a phase of 0. 

% enter spikes with known phases, use a different phase for every trial
Spike = [];
[Spike.trial{1},Spike.time{1}, Spike.timestamp{1}] = deal([]);
rho = linspace(0,1,NTRIALS);
for iTrial = 1:NTRIALS
    SPKFREQ = SINFREQ;
    isi = 1/SPKFREQ; % in that way we will have a perfect phase.
    
    % in cycles
    phase = rho(iTrial);
    nSpikes = 2/isi; % 95 spikes per trials
    spikeSteps = ones(1,nSpikes)*isi;% + 2*rand(1,nSpikes)
    spikeTim = 0 + cumsum(spikeSteps) + phase*isi;
    % find the number of samples
    spikesamples = round(spikeTim/(1/data.fsample))+1; 
    
    spikesamples - (spikeTim/(1/data.fsample)+1) 
   
    spikesamples(spikesamples>length(time{iTrial})) = [];
    trial{iTrial}(3,:) = 0;
    trial{iTrial}(3,spikesamples) = 1;
    
    Spike.time{1}                = [Spike.time{1} spikeTim];
    Spike.trial{1}               = [Spike.trial{1}; iTrial*ones(length(spikeTim),1)];
    Spike.trialtime(iTrial,:)         = [time{1}(1) time{1}(end)];
    Spike.timestamp{1}           = [Spike.timestamp{1} spikeTim*1000000];
end

Spike.label = {'spk1'};
Spike.waveform = {'spk1'};
data.time = time;
data.trial = trial;
data.cfg = cfg;
data.hdr = struct([]);

%%

% test with a simple data structure as the input
cfg = [];
cfg.spikechannel = 3;
cfg.timwin       = [-0.1 0.1];
tic
freq1 = ft_spiketriggeredspectrum(cfg,data);
toc
%
cfg = [];
cfg.method = 'convol';
cfg.foi = [SINFREQ];
cfg.t_ftimwin = 10./cfg.foi;
cfg.taper = 'hanning';
cfg.borderspikes = 'yes';
tic
freq2 = ft_spiketriggeredspectrum(cfg,data);
toc
foi = nearest(freq2.freq,SINFREQ);
%
foi = nearest(freq1.freq,SINFREQ);
inTrial = freq1.trial{1}>0;
phases1 = angle(freq1.fourierspctrm{1}(inTrial,1,foi));
figure, plot(freq1.time{1}(inTrial),phases1,'bx')
%
foi = nearest(freq2.freq,SINFREQ);
inTrial = freq2.trial{1}>0;
phases1 = angle(freq2.fourierspctrm{1}(inTrial,1,foi));
hold on
plot(freq2.time{1}(inTrial),phases1,'ro-')
pause
close 
%%
% test with a simple data structure as the input
cfg = [];
cfg.spikechannel = 3;
cfg.timwin       = [-0.1 0.1];
tic
freq1 = ft_spiketriggeredspectrum(cfg,data);
toc
%
cfg = [];
cfg.method = 'convol';
cfg.foi = [SINFREQ];
cfg.t_ftimwin = 10./cfg.foi;
cfg.taper = 'hanning';
cfg.borderspikes = 'no';
tic
freq2 = ft_spiketriggeredspectrum(cfg,data);
toc
%
cfg = [];
cfg.method = 'convol';
cfg.foi = [SINFREQ];
cfg.t_ftimwin = 10./cfg.foi;
cfg.taper = 'hanning';
cfg.borderspikes = 'no';
cfg.channel = data.label(1:2);
tic
freq3 = ft_spiketriggeredspectrum(cfg,data,Spike);
toc
%%


foi = nearest(freq2.freq,SINFREQ);
%
foi = nearest(freq1.freq,SINFREQ);
inTrial = freq1.trial{1}>0;
phases1 = angle(freq1.fourierspctrm{1}(inTrial,1,foi));
figure, plot(freq1.time{1}(inTrial),phases1,'bx')
%
foi = nearest(freq2.freq,SINFREQ);
inTrial = freq2.trial{1}>0;
phases1 = angle(freq2.fourierspctrm{1}(inTrial,1,foi));
hold on
plot(freq2.time{1}(inTrial),phases1,'ro-')
hold on
foi = nearest(freq3.freq,SINFREQ);
inTrial = freq3.trial{1}>0;
phases1 = angle(freq3.fourierspctrm{1}(inTrial,1,foi));
plot(freq3.time{1}(inTrial),phases1,'gd-')

pause
close 
% NOTE THAT THE PHASE USING THE TRUE TIMES IS MORE RELIABLE!
%%