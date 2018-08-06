function test_bug1262

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_header ft_read_data ft_read_event

% Quick and dirty sanity check for reading of fcdc_buffer_offline data.

% Example data generated with sine example, saved with record.exe:
dirname = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1262/0001');  

hdr = ft_read_header(dirname);
assert(hdr.nChans == 16);
assert(hdr.nSamples == 1216);

dat = ft_read_data(dirname);
assert(all(size(dat) == [hdr.nChans, hdr.nSamples]));

evt = ft_read_event(dirname);
assert(numel(evt) == 0);  % there are no events
