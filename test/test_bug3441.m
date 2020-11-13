function test_bug3441

% MEM 4gb
% WALLTIME 00:20:00
% DEPENDENCY loadcnt

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3441'));

cfg = [];
cfg.datafile = 'mle_pilot.cnt';
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype = 'trigger';
cfg.trialdef.eventvalue =[10];
cfg.trialdef.prestim = .5;
cfg.trialdef.poststim = 1;
cfg = ft_definetrial(cfg);

% STEP 4. PREPROCESS (including baseline correction, re-referencing, and line-noise and band-pass filtering)
cfg.reref='yes'; %turn on re-referencing option
cfg.refchannel={'M1' 'M2'}; %these are the channels to re-reference to
cfg.blc = 'yes'; %include baseline correction
cfg.baselinewindow = [-0.5 0];
cfg.detrend = 'yes';
cfg.dftfilter = 'yes'; %filter out line noise
cfg.dftfreq = 60; %this is the frequency of the line noise (in Europe it would be 50 Hz)
cfg.padding = 10; %amount of seconds to pad each trial to minimize filtering artifacts
cfg.bpfilter='yes'; %bandpass filter
cfg.bpfreq=[.1 100]; %pass frequencies within 0.5 Hz and 90 Hz only
cfg.bpfilttype = 'but'; % or 'fir' or 'firls' (default = 'but')

preproc_odd = ft_preprocessing(cfg);
%save('preproc_odd','preproc_od')

%trial 111, and 113 are the same. both return the same error - "Numeric 
%accuracy issue with the first sample: Rounding off to the nearest integer 
%value". I believe this is coming from loadcnt, see ft bug 2746.
same = isequal(preproc_odd.trial{1,111},preproc_odd.trial{1,113});

% added by Robert: give an error if they are the same
assert(~same);

