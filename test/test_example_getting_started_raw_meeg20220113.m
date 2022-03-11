function test_example_getting_started_with_reading_raw_eeg_or_meg_data

% MEM 4gb
% WALLTIME 00:10:00

%
%% Getting started with reading raw EEG or MEG data
%
% In FieldTrip you first have to define the segments of data in which you are interested, i.e. the "trials". That is done using the DEFINETRIAL function. You can use the DEFINETRIAL function also to show a summary of all events on your data fil
%

datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf');

cfg = [];
cfg.dataset = fullfile(datadir, 'ArtifactMEG.ds');
cfg.trialdef.eventtype  = '?';
ft_definetrial(cfg); % no output variable necessary here, just look at the output in the MATLAB screen

% % The output on screen might look like this
% %
% evaluating trialfunction 'trialfun_general'
% the following events were found in the datafile
% event type: 'trial' with event values:
% no trials have been defined yet, see DEFINETRIAL for further help
% found 76 events
% created 0 trials
% 
% % The important line is
% %
% event type: 'trial' with event values:

% which indicates that the dataset contains 'trials', and that the trials themselves do not have a value. Other events are for example triggers, which usually will have a numeric value.
%
% Since "trial" events have a clear begin, end and duration, you do not have to specify those in DEFINETRIAL. The next step is to call DEFINETRIAL once more, now telling it that it should select the data segments corresponding with the trial events.
%
cfg.trialdef.eventtype = 'trial';
cfg = ft_definetrial(cfg); % now you do want to use an output variable for definetrial, since you need its output

% % This time the following lines will appear on the MATLAB output
% %
% evaluating trialfunction 'trialfun_general'
% found 76 events
% created 76 trials

% which indicate that 76 "fieldtrip-trials" have been made out of the 76 "trial-events" in the datafile. Subsequently you can call the PREPROCESSING function, which will read the desired data segments and (optionally) apply filtering and rereferencing to the data.
%
raw_data = ft_preprocessing(cfg)

% % which will show the following information on screen
% %
% retaining exist trial definition
% retaining exist event information
% found 76 events
% created 76 trials
% rejected    0 trials completely
% rejected    0 trials partially
% resulting  76 trials
% reading and preprocessing
% reading and preprocessing trial 1 from 76
% reading and preprocessing trial 2 from 76
% ...
% reading and preprocessing trial 76 from 76
% 
% raw_data =
%         cfg: [1x1 struct]
%         hdr: [1x1 struct]
%       label: {176x1 cell}
%       trial: {1x76 cell}
%        time: {1x76 cell}
%     fsample: 1200
%        grad: [1x1 struct]

% P.S. prior to calling PREPROCESSING, you might want to do artifact detection using the ARTIFACT_xxx functions (where xxx is for example EOG) and the REJECTARTIFACT function.
%
%% # Another example using "trigger" events
%
% The following example takes the same steps, but in this case the dataset is recorded in pseudo-continuous mode, i.e. there are no gaps between the trials in the data (see \* below). Now we want to read data segments around trigger events.
%
% First determine all event types that are present in the datafile.
%
cfg = [];
cfg.dataset = fullfile(datadir, 'Subject01.ds');
cfg.trialdef.eventtype  = '?';
ft_definetrial(cfg); % no output, just look at screen

% % This will show the following information on scree
% %
% evaluating trialfunction 'trialfun_general'
% the following events were found in the datafile
% event type: 'FC' with event values:
% event type: 'FIC' with event values:
% event type: 'IC' with event values:
% event type: 'STIM' with event values: 16384    65536   196608   327680   589824  1048576
% event type: 'Tr15' with event values:
% event type: 'Tr21' with event values:
% event type: 'Trial' with event values:
% event type: 'backpanel trigger' with event values: 1   3   5   9  16
% event type: 'classification' with event values: 'BAD'
% event type: 'frontpanel trigger' with event values: 16384
% event type: 'trial' with event values:
% no trials have been defined yet, see DEFINETRIAL for further help
% found 1343 events
% created 0 trials

% As you can see, this MEG dataset contains a lot of different events. The interesting events are the "backpanel" triggers. Those triggers can have different values, which are also displayed here. Now specify the trigger value, pre-trigger and post-trigger duration with
%
cfg.trialdef.eventtype  = 'backpanel trigger';
cfg.trialdef.eventvalue = 1;
cfg.trialdef.prestim    = 0.2;
cfg.trialdef.poststim   = 0.6;
cfg = ft_definetrial(cfg); % now remember the output

% % which results in the following info on screen.
% %
% evaluating trialfunction 'trialfun_general'
% found 1343 events
% created 5 trials

% Subsequently, you can do
%
data_raw = ft_preprocessing(cfg)

% % which reads the data and gives following information on screen.
% %
% retaining exist trial definition
% retaining exist event information
% found 1343 events
% created 5 trials
% rejected    0 trials completely
% rejected    0 trials partially
% resulting   5 trials
% reading and preprocessing
% reading and preprocessing trial 1 from 5
% reading and preprocessing trial 2 from 5
% reading and preprocessing trial 3 from 5
% reading and preprocessing trial 4 from 5
% reading and preprocessing trial 5 from 5
% 
% data_raw =
%         cfg: [1x1 struct]
%         hdr: [1x1 struct]
%       label: {187x1 cell}
%       trial: {[187x300 double]  [187x300 double]  [187x300 double]  [187x300 double]  [187x300 double]}
%        time: {[1x300 double]  [1x300 double]  [1x300 double]  [1x300 double]  [1x300 double]}
%     fsample: 300
%        grad: [1x1 struct]
