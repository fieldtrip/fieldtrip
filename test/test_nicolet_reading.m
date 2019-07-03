function test_nicolet_reading

% MEM 2gb
% WALLTIME 00:30:00
% DEPENDENCY test_nicolet_reading_onefile read_nervus_header

% function to test reading of Nicolet/Nervus EEG files
% one is a file from 2006 (older Nicolet format)
% one is a file from 2018 (newer Nicolet format)

path_to_load = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/nicolet');

%%

file1 = 'Patient59_EEG-OPPTAKER-1_t1-NICOLET.e';
file1ascii = 'Patient59_EEG-OPPTAKER-1_t1-ASCII.txt';
test_nicolet_reading_onefile(path_to_load,file1,file1ascii,500,24,45501,1,datetime(2018,06,18,08,07,24));

%%

file2 = 'janbrogger.e';
file2ascii = 'janbrogger.txt';
test_nicolet_reading_onefile(path_to_load,file2,file2ascii,256,24,307201,1,datetime(2006,06,09,13,42,41));

%%
% Test code requested by Robert. For file 1
filepath1 = fullfile(path_to_load,file1);

hdr = ft_read_header(filepath1);
event = ft_read_event(filepath1);
dataopts = {};
alldata = ft_read_data(filepath1, 'header', hdr, dataopts{:});

% open databrowser for file 1
cfg            = [];
cfg.dataset    = filepath1;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data);

%%
% Tests code requested by Robert. For file 2
filepath2 = fullfile(path_to_load,file2);

hdr = ft_read_header(filepath2);
event = ft_read_event(filepath2);
dataopts = {};
alldata = ft_read_data(filepath2, 'header', hdr, dataopts{:});

% open databrowser for file 2
cfg            = [];
cfg.dataset    = filepath2;
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data);
