function test_bug629

% MEM 1gb
% WALLTIME 00:03:01

% TEST test_bug629
% TEST ft_read_header ft_read_data ft_read_event read_mff_header read_mff_data

% note that read_mff_event does not exist, handling of the events is coded in ft_read_event

datadir = '/Volumes/Data/roboos/data/dataformat/testdata/egi/egi_mff Ingrid';

if ~isdir(datadir)
  warning('the test script could not run because the test data is not present');
  return
end

cd(datadir);
dataset = 'pilot05_test 20110120 1433.mff';

hdr   = ft_read_header(dataset);
event = ft_read_event(dataset);
