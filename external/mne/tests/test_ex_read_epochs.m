function test_ex_read_epochs()
% Test IO with FIF epochs

pathstr = fileparts(mfilename('fullpath'));
fname = [pathstr filesep 'data' filesep 'test_raw.fif'];
eventname = [pathstr filesep 'data' filesep 'test-eve.fif'];

event = 1
tmin = -0.2
tmax = 0.5

[data, times, ch_names] = mne_ex_read_epochs(fname, event, eventname, tmin, tmax);

assert(376 == length(ch_names))
assert(7 == length(data))
assert(376 == size(data(1).epoch, 1))
% assert(length(times) == size(data(1).epoch, 2)) # XXX should not fail

