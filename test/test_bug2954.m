function test_bug2954

% WALLTIME 00:10:00
% MEM 2gb

%%
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

%%
clear all
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/edf/shhs1-200001.edf');

%% read the largest subset of channels with a consistent sampling rate, this should return 5 channels at 1Hz
hdr = ft_read_header(filename);
dat = ft_read_data(filename, 'begsample', 1, 'endsample', 100);

assert(numel(hdr.label)==5);
assert(all(size(dat)==[5 100]));

%% read a manually defined subset of channels, this is the set at 1Hz
hdr1  = ft_read_header(filename, 'chanindx', 1);
hdr2  = ft_read_header(filename, 'chanindx', 2);
hdr11 = ft_read_header(filename, 'chanindx', 11);
hdr12 = ft_read_header(filename, 'chanindx', 12);
hdr14 = ft_read_header(filename, 'chanindx', 14);

% read the first channel from each subset
dat1  = ft_read_data(filename, 'header', hdr1 , 'begsample', 1, 'endsample', 100);
dat2  = ft_read_data(filename, 'header', hdr2 , 'begsample', 1, 'endsample', 100);
dat11 = ft_read_data(filename, 'header', hdr11, 'begsample', 1, 'endsample', 100);
dat12 = ft_read_data(filename, 'header', hdr12, 'begsample', 1, 'endsample', 100);
dat14 = ft_read_data(filename, 'header', hdr14, 'begsample', 1, 'endsample', 100);

assert(all(dat1 ==dat(1,:)));
assert(all(dat2 ==dat(2,:)));
assert(all(dat11==dat(3,:)));
assert(all(dat12==dat(4,:)));
assert(all(dat14==dat(5,:)));

assert(isequal(hdr1.label,  hdr.label(1)));
assert(isequal(hdr2.label,  hdr.label(2)));
assert(isequal(hdr11.label, hdr.label(3)));
assert(isequal(hdr12.label, hdr.label(4)));
assert(isequal(hdr14.label, hdr.label(5)));

%% read a manually defined subset of channels, this is the set at 125Hz
hdr125 = ft_read_header(filename, 'chanindx', find(hdr.orig.SampleRate==125));
dat125 = ft_read_data(filename, 'begsample', 1, 'endsample', 1000, 'header', hdr125);

%% the 125Hz channels have the highest sampling rate, the others will be upsampled
data = edf2fieldtrip(filename);

assert(isequal(data.label([3 4 5 8]), hdr125.label));
assert(isequal(data.trial{1}(3,1:1000), dat125(1,:)));
assert(isequal(data.trial{1}(4,1:1000), dat125(2,:)));
assert(isequal(data.trial{1}(5,1:1000), dat125(3,:)));
assert(isequal(data.trial{1}(8,1:1000), dat125(4,:)));



