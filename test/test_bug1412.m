function test_bug1412

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_header loadcnt
% DATA private

hdr = ft_read_header(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1412/Sub1_1.cnt'));

assert(hdr.nSamples == 1505200);

