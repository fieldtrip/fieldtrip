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
trial(1:NTRIALS) = {[cos(SINFREQ*2*pi*time{1});sin(SINFREQ*2*pi*time{1})]}; % so a spike at t = 0 should mean a phase of 0. 

% enter spikes with known phases, use a different phase for every trial
Spike = [];
[Spike.trial{1},Spike.time{1}] = deal([]);
rho = linspace(0,1,20);
for iTrial = 1:NTRIALS
    SPKFREQ = 50;
    isi = 1/SPKFREQ; % in that way we will have a perfect phase.
    
    % in cycles
    phase = rho(iTrial);
    nSpikes = 2/isi - 5; % 95 spikes per trials
    spikeSteps = ones(1,nSpikes)*isi;% + 2*rand(1,nSpikes)
    spikeTim = 0 + cumsum(spikeSteps) + phase*isi;
    % find the number of samples
    spikesamples = round(spikeTim/(1/data.fsample))+1;    
    trial{iTrial}(3,:) = 0;
    trial{iTrial}(3,spikesamples) = 1;
    
    Spike.time{1}                = [Spike.time{1} spikeTim];
    Spike.trial{1}               = [Spike.trial{1}; iTrial*ones(length(spikeTim),1)];
    Spike.trialtime(iTrial,:)         = [time{1}(1) time{1}(end)];
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
tic
freq1 = ft_spiketriggeredspectrum(cfg,data);
foi = nearest(freq1.freq,50);
phases1 = angle(freq1.fourierspctrm{1}(:,1,foi));
figure, plot(phases1)

%%

