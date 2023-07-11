function test_mayo_mef

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY read_mayo_mef21 read_mayo_mef30

%% version 2.1
dataset = '/home/common/matlab/fieldtrip/data/test/pull1378/sample_mef/mef_2p1';

% low-level
hdr = ft_read_header(dataset);
dat = ft_read_data(dataset);
event = ft_read_event(dataset);

% high-level
cfg = [];
cfg.dataset = dataset;
data = ft_preprocessing(cfg);


%% version 3.0
dataset = '/home/common/matlab/fieldtrip/data/test/pull1378/sample_mef/mef_3p0.mef';

% low-level
hdr = ft_read_header(dataset);
dat = ft_read_data(dataset);
event = ft_read_event(dataset);

% high-level
cfg = [];
cfg.dataset = dataset;
data = ft_preprocessing(cfg);
