function test_bug1808

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_header ft_read_sens mne2grad

dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1808/reduced.fif');

hdr  = ft_read_header(dataset);
grad = ft_read_sens(dataset);


