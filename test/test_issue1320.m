function test_issue1320

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY read_edf ft_read_header ft_read_data

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/issue1320/000810000.230120.122109.Signals.Raw.edf');

%%
% this will give a warning, and 27 consistent channels are selected

hdr = ft_read_header(filename);
orig = hdr.orig;

%%
% read all 31 channels

hdr = {};
for i=1:orig.NS
  hdr{i} = ft_read_header(filename, 'chanindx', i);
end

dat = {};
for i=1:orig.NS
  dat{i} = ft_read_data(filename, 'header', hdr{i}, 'chanindx', 1, 'begsample', 1, 'endsample', hdr{i}.Fs);
end

%%
% another way that should work is

data_upsampled = edf2fieldtrip(filename);


