function test_bug629

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_header ft_read_data ft_read_event read_mff_header read_mff_data

% note that read_mff_event does not exist, handling of the events is coded in ft_read_event

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/bug629');

cd(datadir);
dataset = 'pilot05_test 20110120 1433.mff';

hdr   = ft_read_header(dataset);
event = ft_read_event(dataset);
