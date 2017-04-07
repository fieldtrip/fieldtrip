function test_bug1359

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_header ft_read_data loadcnt

% the data is shared with test script 1490
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1490/sub1E3a.cnt');

hdr = ft_read_header(filename);
dat = ft_read_data(filename);

% The report is that the data looks weird and does not contain any
% negative values. I have not checked with old reading functions, but
% with FieldTrip revision 5913 (just after fixing the related bug 1490)
% the data looks meaningful. The channels have large offsets, some appear
% to have clipped, but there are channels that show a physiological signal
% with electrode drift. So it seems that the problem has disappeared.

assert(any(dat(:)>0) && any(dat(:)<0));

