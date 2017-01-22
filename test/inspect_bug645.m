function inspect_bug645

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug645
% TEST ft_read_event

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug645/Events_ECI_TCPIP_55513.xml');

% there is nothing in this script to be tested automatically
% it only serves to determine the timings once

cd /Volumes/Data/roboos/matlab/fieldtrip/fileio/private
tic; ver1 = xml2struct(filename); toc
% xml2struct reading /home/common/matlab/fieldtrip/data/test/bug645/Events_ECI_TCPIP_55513.xml
% Elapsed time is 22.311995 seconds.

cd /Volumes/Data/roboos/matlab/fieldtrip/external/xml4mat
tic; ver2 = xml2struct(filename); toc
% Elapsed time is 404.807848 seconds.


