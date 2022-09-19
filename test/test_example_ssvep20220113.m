function test_example_ssvep

% MEM 4gb
% WALLTIME 00:10:00

%
%% Analyze Steady-State Visual Evoked Potentials (SSVEPs)
%
% Steady-state stimulation is frequently used for sensory stimulation in the visual (SSVEP), auditory (SSAEP), and somatosensory (SSSEP) domains. On this page, we will first present an example analysis strategy for a 64-channel SSVEP dataset. Subsequently, we will make some simulated data and use another analysis strategy.
%
%% # Analyze experimental data
%
% The design of a trial is shown in the figure below. On each trial, a ring is flashing at 8.33 Hz for 4.56 s (at 100 Hz, on for 6 frames and off for 6 frames). During this time, a 100-ms rectangle is shown about every 1 s for four times. In a block of 40 trials, the task is either low load (i.e., count color, e.g., black) or high load (i.e., count a certain combination, e.g., black horizontal rectangle and white vertical rectangle). After each trial, subjects report if they counted 2 or 3 targets. Within each block, the ring is shown at four excentricities: 2, 3, 4, and 6 visual degrees (10 trials at each excentricity), in pseudorandom order for each subject. Load levels alternate between blocks (i.e., LHLHLHLH or HLHLHLHL , counterbalanced across subjects).
%
% There are 2 loads x 4 blocks x 4 ring excentricities x 10 trials = 320 SSVEP trials. Because 4 rectangles are shown within each trial, this gives 4 x 320 = 1280 rectangle trials. This gives a total of 320 + 1280 = 1600 trials.
%
%
% In the analysis, consider whether the frequency tagged (or steady state ) stimulus is
%
% #  phase consistent over trials. If so, average and then do wavelet/mtmconvol/mtmfft, or time-domain regression. (In the example, each 4.56-s SSVEP trial has the same phase.)
% #  not phase consistent across trials, but the phase of the stimulus is known. If so, Fourier decompose to get the complex representation, deal with single-trial phase differences, then average. (In the example, the onset of the rectangles is jittered relative to the onset of a SSVEP trial. Thus, the phase of the flashing rings varies between rectangle trials).
% #  not phase consistent across trials, and the phase of the stimulus is not known. If so, Fourier decompose to get the power and average over trials.
%
% Furthermore, consider whether the cortical response
%
% #  is assumed to be constant within the trial
% #  changes over time within the trial
%
% Finally, consider whether the stimulation contains
%
% #  a single frequency, as in a traditional SSVEP
% #  a mixture of multiple frequencies, as in frequency tagging
%
%% ## Using time-domain analysis
%
% This section is still to be written.
%
%% ## Using frequency analysis
%
% This section should be made specific to the example dataset.
%
%% # Create and analyze a simulated steady-state dataset
%
fsample = 1000;
nsample = 100*fsample;

% start with a continuous representation of 100 seconds of data
data.label = {'trigger', 'eeg'};
data.time = {(1:nsample)/fsample};
data.trial = {zeros(length(data.label), nsample)};

% add a trigger every 100ms
i = 100;
while i<nsample
  data.trial{1}(1,i) = 1;
  i = i + 100;
end

% create some sort of SSVEP signal
data.trial{1}(2,:) = ft_preproc_bandpassfilter(data.trial{1}(1,:), fsample, [3 18], [], [], 'onepass');
data.trial{1}(2,:) = data.trial{1}(2,:) + 0.02*randn(size(data.trial{1}(2,:)));

% cut into one-second snippets
cfg = [];
cfg.length = 1;
data = ft_redefinetrial(cfg, data);

plot(data.time{1}, data.trial{1})
legend(data.label)

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'powandcsd';
cfg.foilim = [1 100];
freq = ft_freqanalysis(cfg, data);

plot(freq.freq, freq.powspctrm);
legend(freq.label)

% normalize the CSD for the power in the trigger sequence
freq.crsspctrm = freq.crsspctrm ./ freq.powspctrm(1,:);

plot(freq.freq, abs(freq.crsspctrm));
xlabel('frequency (Hz)')
ylabel('phase-locked amplitude (a.u.)')
