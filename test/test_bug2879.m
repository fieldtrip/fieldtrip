function test_bug2879

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourcestatistics
% DATA private

filename = dccnpath('/project/3031000.02/test/bug2879.mat');
load(filename)

stat = ft_sourcestatistics(cfg,dat1{:},dat2{:});

