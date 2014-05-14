function test_bug2096

% MEM 2000mb
% WALLTIME 00:10:00

% TEST test_bug2096
% TEST ft_sourcewrite

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2096.mat'));

cfg = [];
cfg.filename = [tempname '.cii'];
ft_sourcewrite(cfg, source);

