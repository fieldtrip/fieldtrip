function test_pull1419

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY


filename = dccnpath('/home/common/matlab/fieldtrip/data/test/pull1419/ANM255201_20141125_trimmed.nwb');

% low-level reading functions
hdr = ft_read_header(filename);
dat = ft_read_data(filename);
event = ft_read_event(filename); % returns nothing at the moment

% high-level reading function
cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);
