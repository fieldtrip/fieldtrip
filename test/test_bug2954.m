% function test_bug2954

% WALLTIME 00:10:00
% MEM 1gb

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/edf/0601_s.edf');

% all channels
hdr = ft_read_header(filename);
dat = ft_read_data(filename, 'begsample', 1, 'endsample', 100);


%% read the first channel from the header
hdr1 = ft_read_header(filename, 'chanindx', 1);
hdr2 = ft_read_header(filename, 'chanindx', 2);
hdr3 = ft_read_header(filename, 'chanindx', 3);
hdr4 = ft_read_header(filename, 'chanindx', 4);
dat1 = ft_read_data(filename, 'chanindx', 1, 'header', hdr1, 'begsample', 1, 'endsample', 100);
dat2 = ft_read_data(filename, 'chanindx', 1, 'header', hdr2, 'begsample', 1, 'endsample', 100);
dat3 = ft_read_data(filename, 'chanindx', 1, 'header', hdr3, 'begsample', 1, 'endsample', 100);
dat4 = ft_read_data(filename, 'chanindx', 1, 'header', hdr4, 'begsample', 1, 'endsample', 100);

assert(all(dat1==dat(1,:)));
assert(all(dat2==dat(2,:)));
assert(all(dat3==dat(3,:)));
assert(all(dat4==dat(4,:)));

%% read all channels from the header, i.e. one channel each
hdr1 = ft_read_header(filename, 'chanindx', 1);
hdr2 = ft_read_header(filename, 'chanindx', 2);
hdr3 = ft_read_header(filename, 'chanindx', 3);
hdr4 = ft_read_header(filename, 'chanindx', 4);
dat1 = ft_read_data(filename, 'header', hdr1, 'begsample', 1, 'endsample', 100);
dat2 = ft_read_data(filename, 'header', hdr2, 'begsample', 1, 'endsample', 100);
dat3 = ft_read_data(filename, 'header', hdr3, 'begsample', 1, 'endsample', 100);
dat4 = ft_read_data(filename, 'header', hdr4, 'begsample', 1, 'endsample', 100);

assert(all(dat1==dat(1,:)));
assert(all(dat2==dat(2,:)));
assert(all(dat3==dat(3,:)));
assert(all(dat4==dat(4,:)));

%% read the selected channels from the full header
hdr = ft_read_header(filename);
dat1 = ft_read_data(filename, 'chanindx', 1, 'header', hdr, 'begsample', 1, 'endsample', 100);
dat2 = ft_read_data(filename, 'chanindx', 2, 'header', hdr, 'begsample', 1, 'endsample', 100);
dat3 = ft_read_data(filename, 'chanindx', 3, 'header', hdr, 'begsample', 1, 'endsample', 100);
dat4 = ft_read_data(filename, 'chanindx', 4, 'header', hdr, 'begsample', 1, 'endsample', 100);

assert(all(dat1==dat(1,:)));
assert(all(dat2==dat(2,:)));
assert(all(dat3==dat(3,:)));
assert(all(dat4==dat(4,:)));

