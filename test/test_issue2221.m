function test_issue2221

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY read_biosemi_bdf
% DATA private

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/issue2221.bdf');

hdr = ft_read_header(filename);
dat = ft_read_data(filename);
evt = ft_read_event(filename);

% there are 66 blocks of 1 second in the file
assert(hdr.nSamples/4096==66);
