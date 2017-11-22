function test_bug2887

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_read_header ft_read_data ft_read_event read_edf

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2887/EDFtest.edf');

hdr = ft_read_header(filename);
dat = ft_read_data  (filename);
evt = ft_read_event (filename);


