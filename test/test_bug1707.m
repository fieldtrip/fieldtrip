function test_bug1707

% MEM 2gb
% WALLTIME 0:05:00

% TEST test_bug1707
% TEST ft_channelrepair ft_combineplanar ft_componentanalysis
% TEST ft_connectivityanalysis ft_freqanalysis ft_freqdescriptives
% TEST ft_freqgrandaverage ft_freqstatistics ft_megplanar
% TEST ft_prepare_neighbours ft_preprocessing ft_singleplotER
% TEST ft_singleplotTFR ft_sourceanalysis ft_sourcedescriptives
% TEST ft_sourcegrandaverage ft_sourceplot ft_timelockanalysis
% TEST ft_timelockgrandaverage ft_topoplotER ft_topoplotCC ft_topoplotIC
% TEST ft_multiplotCC ft_multiplotER ft_multiplotTFR

% bug description:
% Anecdotally, Stephen reported issues relating to different channels and
% different orderings of channels that determine the output of
% ft_timelockstatistics, depending on what the order of the input arguments has
% been.
%
% FieldTrip should nowhere assume that the order and list of channels present is
% the same for all data input arguments, nor that the order is the same in e.g.
% the data.label and data.grad.label.
% 
% TODO; make a test script that inputs stuff and outputs stuff
%
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1707


% all data test cases
datainfo = ref_datasets;
%sel      = match_str({datainfo.datatype}',{'ctf275' 'itab153' 'yokogawa160'}');
%datainfo = datainfo(sel);

% all functions + appropriate cfgs to be tested

%% check whether all runs through smoothly

if false
% I (=robert) noticed that the code below is still rather incomplete
% rather than it showing up as error on the dashboard, I have now commented it out for the time being

% first step: preprocess all the data
% ft_channelrepair
sel      = match_str({datainfo.datatype}',{'bdf' 'brainvision' 'edf'}');
data     = datainfo(sel)
funs{end+1} = 'ft_channelrepair';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_combineplanar';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_componentanalysis';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_connectivityanalysis';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_freqanalysis';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_freqdescriptives';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_freqgrandaverage';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_freqstatistics';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_megplanar';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_prepare_neighbours';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_preprocessing';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_singleplotER';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_singleplotTFR';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_sourceanalysis';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_sourcedescriptives';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_sourcegrandaverage';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_sourceplot';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_timelockanalysis';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_timelockgrandaverage';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_topoplotER';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_topoplotCC';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_topoplotIC';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_multiplotCC';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_multiplotER';
cfgs{end+1} = [];
cfgs{end}.

funs{end+1} = 'ft_multiplotTFR';
cfgs{end+1} = [];
cfgs{end}.

end
