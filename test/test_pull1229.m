function test_pull1229

% MEM 8gb
% WALLTIME 00:20:00
% DEPENDENCY ft_read_header ft_read_data ft_read_event

% this is dataset "d-18-s-3" from https://osf.io/52gy7/
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/pull1229'));

d = dir('*.edf');
filename = {d.name};

hdr = ft_read_header(filename);

dat = ft_read_data(filename);
dat = ft_read_data(filename, 'header', hdr);
dat = ft_read_data(filename, 'header', hdr, 'begsample', 1, 'endsample', 2560);

event = ft_read_event(filename);
event = ft_read_event(filename, 'header', hdr);
