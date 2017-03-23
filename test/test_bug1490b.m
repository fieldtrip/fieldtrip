function test_bug1490b

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_read_header ft_read_data loadcnt

% this is a second test pertaining to http://bugzilla.fcdonders.nl/show_bug.cgi?id=1490#c11

warning('this bug cannot be fixed in FieldTrip, but should be fixed in EEGLAB');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1490/cba1ff01.cnt');
dataformat = 'ns_cnt16';

hdr = ft_read_header(filename);

% try to titrate the option related to the error, the first ones all used to fail
dat = ft_read_data(filename, 'header',  hdr, 'dataformat', dataformat, 'begsample',       1, 'endsample', 206440, 'chanindx', 1:hdr.nChans, 'checkboundary', 1, 'dataformat', dataformat);
dat = ft_read_data(filename, 'header',  hdr, 'dataformat', dataformat, 'begsample',  105704, 'endsample', 206440, 'chanindx', 1:hdr.nChans, 'checkboundary', 1, 'dataformat', dataformat);
dat = ft_read_data(filename, 'header',  hdr, 'dataformat', dataformat, 'begsample',  105704, 'endsample', 206440, 'chanindx', 1:hdr.nChans);
dat = ft_read_data(filename, 'header',  hdr, 'dataformat', dataformat, 'begsample',  105704, 'endsample', 206440, 'chanindx', 1);
dat = ft_read_data(filename,                 'dataformat', dataformat, 'begsample',  105704, 'endsample', 206440, 'chanindx', 1);
dat = ft_read_data(filename,                 'dataformat', dataformat, 'begsample',  105704, 'endsample', inf,    'chanindx', 1); % this failed on 22 aug 2012
dat = ft_read_data(filename,                 'dataformat', dataformat, 'begsample',  105704, 'endsample', inf);                   % this failed on 22 aug 2012
dat = ft_read_data(filename,                 'dataformat', dataformat, 'begsample',       1, 'endsample', inf,    'chanindx', 1); % this worked on 22 aug 2012
dat = ft_read_data(filename,                 'dataformat', dataformat, 'begsample',       1, 'endsample', inf);                   % this worked on 22 aug 2012

% this particular file seems to be in 40-sample blocks
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 1, 'endsample', 80);
assert(size(dat,2)==80, 'incorrect number of samples');
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 1, 'endsample', 40);
assert(size(dat,2)==40, 'incorrect number of samples');
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 1, 'endsample', 39);
assert(size(dat,2)==39, 'incorrect number of samples');
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 1, 'endsample', 41);
assert(size(dat,2)==41, 'incorrect number of samples');
dat = ft_read_data(filename, 'dataformat', dataformat, 'begsample', 2, 'endsample', 41);
assert(size(dat,2)==40, 'incorrect number of samples');

