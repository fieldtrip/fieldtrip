function test_bug2170

% WALLTIME 00:03:01

% TEST test_bug2170
% TEST ft_filetype ft_read_event

filename1 = dccnfilename('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/raw.fif');
filename2 = dccnfilename('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/run_01_raw.fif');

assert(ft_filetype(filename1, 'neuromag_fif'));
assert(ft_filetype(filename2, 'neuromag_fif'));


