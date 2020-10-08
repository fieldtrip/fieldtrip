function test_issue1554()

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_spike_maketrials

%% Raw spike structure with 10 spikes randomly distributed in each of a 100 trials (1000 total spikes) 
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
spike.label{1}      = 'spike1';
spike.waveform{1}   = rand(32,length(tsAll));

%% Add amplitude field with same dimensions as timestamps
spike.amplitude{1} = rand(size(tsAll)); % in my case these are floats

% now, for the sake of comparison later, let's have timepoint 15 in every
% waveform take the value that is in amplitude
spike.waveform{1}(15,:) = spike.amplitude{1} ;

%% Make trial structure with all trials where the problem does not occur
cfgC = [];
cfgC.trl = uint64(round(events));
cfgC.timestampspersecond = 1000000;
cfgC.trl(:,3) = -1*1000000;

% spiketime is now organized in trials (66 trials x 10 spikes = 660 spikes)
S1 = ft_spike_maketrials(cfgC, spike);

% Pick a random spike from the timelocked structure
spikenr = randi(size(S1.timestamp{1}, 2));
fprintf('Spike %d is at %0.0f milliseconds in trial %d\n', spikenr, S1.time{1}(spikenr)*1000, S1.trial{1}(spikenr));

% so is the waveform for that matter
figure; 
hold on
title(fprintf('Waveform of spike %d\n', spikenr));
plot(squeeze(S1.waveform{1}(1, :, spikenr)))

% and what about amplitude?
fprintf('Spike %d has an amplitude of %f\n', spikenr, S1.amplitude{1}(spikenr));
fprintf('Spike %d has an amplitude of: %f, the same as its waveform: %f\n', spikenr, S1.amplitude{1}(spikenr), S1.waveform{1}(1, 15, spikenr));


% Now remove some trials so that only a subset of spikes will end up in the
% timelocked structure
cfgC.trl = cfgC.trl(33:100, :);
S2 = ft_spike_maketrials(cfgC, spike);

% Pick a random spike from the new timelocked structure
spikenr = randi(size(S2.timestamp{1}, 2));

% The amplitude field is not adjusted to the new trial structure, so the
% it cannot be indexed in the same way as the time, waveform and timestamp
% field:
if S2.amplitude{1}(spikenr) ~= S2.waveform{1}(1, 15, spikenr)
    error(sprintf('Spike %d has an amplitude of: %f, NOT the same as its waveform: %f\n', spikenr, S2.amplitude{1}(spikenr), S2.waveform{1}(1, 15, spikenr)))
end

% Or another way to put is, is that this works:
[S2.trial{1}; S2.timestamp{1}];

% But this doesn't:
[S2.trial{1}; S2.amplitude{1}];
