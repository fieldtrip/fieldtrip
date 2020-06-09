function test_issue1425

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preprocessing ft_read_data

d = dir(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1425/*.edf'));

hdr = ft_read_header(fullfile(d.folder,d.name));
event = ft_read_event(fullfile(d.folder,d.name));

