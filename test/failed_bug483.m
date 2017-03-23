function failed_bug483

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug483
% TEST ft_read_header ft_read_data ft_read_event

% the original file location was /Users/robert/Manzana/data/dataformat/testdata/neuromag/rik_henson_MRC-CBU/raw.fif
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/raw.fif');

% open a new file and close it again
fidold = fopen(tempname, 'wb');
fclose(fidold);

hdr = ft_read_header(filename);
hdr = ft_read_header(filename);
hdr = ft_read_header(filename);

dat = ft_read_data(filename, 'begsample', 1, 'endsample', 1000, 'ckeckbondary', false);
dat = ft_read_data(filename, 'begsample', 1, 'endsample', 1000, 'ckeckbondary', false);
dat = ft_read_data(filename, 'begsample', 1, 'endsample', 1000, 'ckeckbondary', false);

event = ft_read_event(filename);
event = ft_read_event(filename);
event = ft_read_event(filename);

% open a new file and close it again
fidnew = fopen(tempname, 'wb');
fclose(fidnew);

% the file identifier for a new file should not have changed
% if it changed, then probably there is still a fif file open
assert(fidold==fidnew);


