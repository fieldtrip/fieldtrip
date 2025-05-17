function test_nwb

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_header
% DATA private

% See https://github.com/fieldtrip/fieldtrip/pull/1419
% See https://github.com/fieldtrip/fieldtrip/issues/2499

filename1 = dccnpath('/project/3031000.02/test/original/eeg/nwb/test.nwb');
filename2 = dccnpath('/project/3031000.02/test/original/eeg/nwb/sub-002_task-FaceRecognition_eeg.nwb');

%%

hdr = ft_read_header(filename1);
dat = ft_read_data(filename1);
evt = ft_read_event(filename1);

%%

hdr = ft_read_header(filename2);
dat = ft_read_data(filename2);
evt = ft_read_event(filename2);
