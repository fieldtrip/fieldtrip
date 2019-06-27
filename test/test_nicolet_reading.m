function test_nicolet_reading

% MEM 100mb
% WALLTIME 00:02:00

% function to test reading of Nicolet/Nervus EEG files
% one is a file from 2006 (older Nicolet format)
% one is a file from 2018 (newer Nicolet format)

%path_to_load = dccnpath('/home/common/matlab/fieldtrip');
path_to_load = 'C:\Midlertidig_Lagring\fieldtrip\test\private\test_nicolet_reading\';

test_nicolet_reading_onefile(path_to_load,'Patient59_EEG-OPPTAKER-1_t1-NICOLET.e','Patient59_EEG-OPPTAKER-1_t1-ASCII.txt',500,24,45501,1,datetime(2018,06,18,08,07,24));
test_nicolet_reading_onefile(path_to_load,'janbrogger.e','janbrogger.txt',256,24,307201,1,datetime(2006,06,09,13,42,41));

    
end