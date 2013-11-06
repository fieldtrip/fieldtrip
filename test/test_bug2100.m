function test_bug2100

% MEM 1500mb
% WALLTIME 00:03:06

% TEST test_bug2100
% TEST ft_read_mri read_ctf_mri4

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2100/Sub02.mri');
mri = ft_read_mri(filename);

