function test_bug1348

% MEM 2gb
% WALLTIME 00:10:00

% TEST test_bug1348
% TEST ft_read_header ft_read_data loadcnt

filename = '/home/common/matlab/fieldtrip/data/test/bug1348/E3a.cnt';

hdr = ft_read_header(filename);
disp(hdr);

