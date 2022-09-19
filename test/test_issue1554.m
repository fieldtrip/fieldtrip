function test_issue1554

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_spike_maketrials

clear all 

%% start with a standard raw spike structure with 10 spikes randomly distributed in each of a 100 trials (1000 total spikes)
spike   = [];
tsAll   = [];
events  = [];

for iTrial = 1:100
  timebeg           = 0;
  ts                = sort(uint64(timebeg+iTrial*1000000+round(1000000*rand(1,10))));
  tsAll             = [tsAll; ts(:)];
  events(iTrial,:)  = [timebeg+iTrial*1000000 timebeg+(iTrial+1)*1000000];
end

spike.timestamp{1}  = tsAll;
spike.unit{1}       = ones(size(spike.timestamp{1}));
spike.timestamp{2}  = tsAll;
spike.unit{2}       = ones(size(spike.timestamp{1}));

spike.label         = {'spike1', 'spike2'};
spike.waveform{1}   = rand(1, 32,length(tsAll));
spike.waveform{2}   = rand(1, 32,length(tsAll));

% ensure it is a correct structure
spike = ft_checkdata(spike, 'datatype', 'spike');

%% Add amplitude field with same dimensions as timestamps
spike.amplitude{1} = rand(size(spike.timestamp{1})); % in my case these are floats
spike.amplitude{2} = rand(size(spike.timestamp{2})); % in my case these are floats
spike.amplitudedimord = '{chan}_spike';

% now, for the sake of comparison later, let's have timepoint 15 in every
% waveform take the same value that is in amplitude
spike.waveform{1}(1, 15,:) = spike.amplitude{1};

%% Make trial structure with all trials where the problem does not occur

cfg = [];
cfg.trl = uint64(round(events));
cfg.trl(:,3) = -1*1000000;
cfg.timestampspersecond = 1000000;

% spiketime is now organized in trials (100 trials x 10 spikes = 1000 spikes)
spike_trl = ft_spike_maketrials(cfg, spike);

figure
plot(spike_trl.timestamp{1}, spike_trl.trial{1}, '.'); xlabel('timestamp'); ylabel('trial')

% ensure that the timestamps, waveforms and amplitudes are the same
assert(isequal(spike.timestamp, spike_trl.timestamp));
assert(isequal(spike.amplitude, spike_trl.amplitude));
assert(isequal(spike.waveform,  spike_trl.waveform));

%% Pick a random spike from the timelocked structure
spikenr = randi(size(spike_trl.timestamp{1}, 2));
fprintf('Spike %d is at %0.0f milliseconds in trial %d\n', spikenr, spike_trl.time{1}(spikenr)*1000, spike_trl.trial{1}(spikenr));

% so is the waveform for that matter
figure;
hold on
title(fprintf('Waveform of spike %d\n', spikenr));
plot(squeeze(spike_trl.waveform{1}(1, :, spikenr)))

% and what about amplitude?
fprintf('Spike %d has an amplitude of %f\n', spikenr, spike_trl.amplitude{1}(spikenr));
fprintf('Spike %d has an amplitude of: %f, the same as its waveform: %f\n', spikenr, spike_trl.amplitude{1}(spikenr), spike_trl.waveform{1}(1, 15, spikenr));

%% Now remove some trials so that only a subset of spikes will end up in the timelocked structure
cfg = [];
cfg.trl = uint64(round(events));
cfg.trl(:,3) = -1*1000000;
cfg.timestampspersecond = 1000000;
cfg.trl = cfg.trl(33:66, :);
spike_trl = ft_spike_maketrials(cfg, spike);

%% Pick a random spike from the timelocked structure
spikenr = randi(size(spike_trl.timestamp{1}, 2));
fprintf('Spike %d is at %0.0f milliseconds in trial %d\n', spikenr, spike_trl.time{1}(spikenr)*1000, spike_trl.trial{1}(spikenr));

% so is the waveform for that matter
figure;
hold on
title(fprintf('Waveform of spike %d\n', spikenr));
plot(squeeze(spike_trl.waveform{1}(1, :, spikenr)))

% and what about amplitude?
fprintf('Spike %d has an amplitude of %f\n', spikenr, spike_trl.amplitude{1}(spikenr));
fprintf('Spike %d has an amplitude of: %f, the same as its waveform: %f\n', spikenr, spike_trl.amplitude{1}(spikenr), spike_trl.waveform{1}(1, 15, spikenr));

