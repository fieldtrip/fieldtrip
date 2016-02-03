function test_bug1207

% MEM 1gb
% WALLTIME 00:10:00

% TT5.ntt in this directory
% TT5_1.t in this directory

spike1 = ft_read_spike('TT5.ntt');
spike2 = ft_read_spike('TT5_1.t');
