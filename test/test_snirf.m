function test_snirf

% WALLTIME 00:10:00
% MEM 2GB
% DEPENDENCY homer2fieldtrip fieldtrip2homer ft_write_data

p = tempdir;
f1 = fullfile(p, 'data1.snirf');
f2 = fullfile(p, 'data2.snirf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test the conversion of a Homer file to snirf

% the Homer nirs file has the 's' channel, which is not present in the snirf representation

% read the original data
cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/original/nirs/homer/S1001_run01.nirs');
cfg.channel = {'all', '-s*'};
data1 = ft_preprocessing(cfg);
event1 = ft_read_event(cfg.dataset);

% write the data to disk and convert to snirf
chanindx = find(ismember(data1.hdr.label, data1.label));
ft_write_data(f1, data1.trial{1}, 'header', data1.hdr, 'chanindx', chanindx, 'event', event1)

% read the converted data
cfg = [];
cfg.dataset = f1;
data1a = ft_preprocessing(cfg);
event1a = ft_read_event(cfg.dataset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test the conversion of an Artinis file to snirf

% read the original data
cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/original/nirs/artinis/Helena/190528_fingertap_L.oxy3');
cd(fileparts(cfg.dataset)); % we should be located in the directory with the optodetemplates.xml
data2 = ft_preprocessing(cfg);
event2 = ft_read_event(cfg.dataset, 'chanindx', -1); % do not parse the ADC channels

% replace the original data of the first sample by channel numbers. This helps to assert
% whether data2a is still the same as the original data (data2), even if the
% channel order is switched.
data2.trial{1}(:,1) = [1:length(data2.label)]';

% write the data to disk and convert to snirf
ft_write_data(f2, data2.trial{1}, 'header', data2.hdr, 'event', event2)

% read the converted data
cfg = [];
cfg.dataset = f2;
data2a = ft_preprocessing(cfg);
event2a = ft_read_event(cfg.dataset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up

delete(f1)
delete(f2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare all fields

assert(isalmostequal(data1.time,        data1a.time, 'abstol', 1e-8));
assert(isalmostequal(data1.fsample,     data1a.fsample, 'abstol', 1e-8));
assert(isalmostequal(data1.sampleinfo,  data1a.sampleinfo, 'abstol', 1e-8));
assert(isequal(data1.label, data1a.label));
assert(isalmostequal(data1.trial,       data1a.trial, 'abstol', 1e-8));

% the duration of the Homer and snirf events is not the same
assert(isequal({event1.type},   {event1a.type}));
assert(isequal([event1.sample], [event1a.sample]));
assert(isequal([event1.value],  [event1a.value]));

%%

% do not compare the label field for Artinis, it differs in Tx/Rx versus S/D
assert(isalmostequal(data2.time,        data2a.time, 'abstol', 1e-8));
assert(isalmostequal(data2.fsample,     data2a.fsample, 'abstol', 1e-8));
assert(isalmostequal(data2.sampleinfo,  data2a.sampleinfo, 'abstol', 1e-8));
assert(length(intersect(data2.trial{1}(:,1), data2a.trial{1}(:,1)))==length(data2.trial{1}(:,1))); % the first sample should contain 1:nchan (see here above)

% the channel names are not the same, neither for the nirs channels, not for the ADC/AUX channels
% assert(isequal(data2.label, data2a.label));

% the type and value of the Artinis and snirf events are not the same
assert(isequal([event2.sample], [event2a.sample]));
