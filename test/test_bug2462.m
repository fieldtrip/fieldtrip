function test_bug2462

% MEM 1500mb
% WALLTIME 00:10:00
% test_bug2462
cfg = [];
cfg.dataset = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2462/scan1_Filters_125HzLP.dat');

data = ft_preprocessing(cfg);
end


