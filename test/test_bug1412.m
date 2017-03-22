function test_bug1412

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_read_header loadcnt

hdr = ft_read_header(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1412/Sub1_1.cnt'));

assert(hdr.nSamples == 1505200);

