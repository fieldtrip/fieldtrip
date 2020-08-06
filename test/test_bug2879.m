function test_bug2879

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourcestatistics

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2879');
load(filename)

stat = ft_sourcestatistics(cfg,dat1{:},dat2{:});

