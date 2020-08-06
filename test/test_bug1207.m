function test_bug1207

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY

spike1 = ft_read_spike(dccnpath('/home/common/matlab/fieldtrip/data/test/spike/TT5.ntt'));
spike2 = ft_read_spike(dccnpath('/home/common/matlab/fieldtrip/data/test/spike/TT5_1.t'));
