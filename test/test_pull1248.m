function test_pull1248

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_header ft_read_data ft_read_event dataset2files

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/pull1248/1106miedo.vhdr');

%%

hdr = ft_read_header(filename);
dat = ft_read_data(filename);
evt = ft_read_event(filename); % this used to crash due to there not being a vmrk file
