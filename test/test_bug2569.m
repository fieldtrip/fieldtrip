function test_bug2569

% WALLTIME 00:10:00
% MEM 6gb

% TEST read_wdq_header read_wdq_data read_wdq_event

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2569/WoD dec2013 rat nr 11.WDQ');
ft_databrowser(cfg)

hdr = ft_read_header(cfg.dataset);
dat = ft_read_data(cfg.dataset);
evt = ft_read_event(cfg.dataset);

assert(hdr.nChans==8);
assert(size(dat,1)==8);
assert(size(dat,2)==20408500);
assert(length(evt)==12);
