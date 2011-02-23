function test_ex_read_raw()
% Test IO with FIF raw files

pathstr = fileparts(mfilename('fullpath'));
fname = [pathstr filesep 'data' filesep 'test_raw.fif'];

from = 50;
to = 55;
in_samples = false;
[data, times] = mne_ex_read_raw(fname, from, to, in_samples);

assert(size(data, 2) == length(times))
assertTrue(all(times < 56))
assertTrue(all(times > 49))

raw = fiff_setup_read_raw(fname);
fiff_write_raw_segment_times('foo.fif', raw, 51, 53);

raw = fiff_setup_read_raw('foo.fif');
[data, times] = fiff_read_raw_segment(raw);
assert(size(data, 2) == length(times))
assertTrue(all(times < 53.01))
assertTrue(all(times > 50.99))

