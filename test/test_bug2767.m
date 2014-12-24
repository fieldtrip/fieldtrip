function test_bug2767

% WALLTIME 00:10:00
% MEM 1gb

% TEST test_bug2767

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/01_ljh_firststd_meg_182.fif');

hdr = ft_read_header(filename);
dat = ft_read_data(filename);
evt = ft_read_event(filename);

cfg = [];
cfg.dataset = filename;
raw = ft_preprocessing(cfg);

cfg = [];
erf = ft_timelockanalysis(cfg, raw);

