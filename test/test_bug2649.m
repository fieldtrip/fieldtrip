function test_bug2649

% WALLTIME 00:10:00
% MEM 2gb

% TEST ft_write_data write_brainvision_eeg

fileorig = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/brainvision/Mischa.eeg');

orig_hdr = ft_read_header(fileorig);
orig_dat = ft_read_data(fileorig);
orig_evt = ft_read_event(fileorig);

filecopy = [tempname '.eeg'];

ft_write_data(filecopy, orig_dat, 'header', orig_hdr, 'event', orig_evt, 'dataformat', 'brainvision_eeg');

copy_hdr = ft_read_header(filecopy);
copy_dat = ft_read_data(filecopy);
copy_evt = ft_read_event(filecopy);

assert(isalmostequal(orig_dat, copy_dat, 'reltol', 1e-6));
assert(isequal(keepfields(orig_hdr, {'Fs', 'nChans', 'label', 'nSamples', 'nSampelsPre'}), keepfields(copy_hdr, {'Fs', 'nChans', 'label', 'nSamples', 'nSampelsPre'})));
assert(isequal(keepfields(orig_evt, {'type', 'value', 'sample', 'duration'}), keepfields(copy_evt, {'type', 'value', 'sample', 'duration'})));

% clean up the temp files
[p, f, x] = fileparts(filecopy);
delete(fullfile(p, [f '.eeg']));
delete(fullfile(p, [f '.vhdr']));
delete(fullfile(p, [f '.vmrk']));
