function test_bug2100

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_mri read_ctf_mri4

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2100/Sub02.mri');
mri = ft_read_mri(filename);

