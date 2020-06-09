%% Lag extraction tutorial script
% this script runs an EEGLAB plugin used for lag extraction on single trial data.
%
% $Id: lagextraction_tutorial.m 4 2009-08-15 21:10:35Z gramfort $

%% Load data
load('data/oddball3-num1-512Hz-chan10.set','-mat');

%% Set parameters and run lag extraction
use_ica = false; % Set to true, if you want to realign based on an ICA component
channel = 1; % Index of channel or ICA component used for realignment
time_win = [150 500]; % (ms) : work on this time window
bad_trials = []; % set bad trials

clear options
options.sigma = [0.01:0.01:0.2];
options.alpha = [0.001,0.01,0.1];
options.disp_log = false;
[EEG, com, order, lags, event_type, E_lags] = pop_extractlag( EEG , use_ica, channel, time_win, options);

%% View ERP image reordered
figure;
pop_erpimage(EEG,1, [channel],[],EEG.chanlocs(channel).labels,1,1,{ event_type }, ...
             [],'latency' ,'yerplabel','\muV','erp','cbar');

%% Re-epoch the data
EEG = pop_epoch( EEG, {event_type}, [-0.4 0.3]);

%% View ERP image of re-epoched data
figure;
pop_erpimage(EEG,1, [channel],[],EEG.chanlocs(channel).labels,1,1,{},[],'','yerplabel','\muV','erp','cbar');