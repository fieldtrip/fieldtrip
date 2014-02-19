function test_bug1412

% MEM 2gb
% WALLTIME 00:10:00

% TEST test_bug1412
% TEST ft_read_header loadcnt

filename = '/home/common/matlab/fieldtrip/data/test/bug1412/Sub1_1.cnt';

hdr = ft_read_header(filename)

assert(hdr.nSamples == 1505200);

