function test_bug1348

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_read_header ft_read_data loadcnt

hdr = ft_read_header(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1348/E3a.cnt'));
disp(hdr);

