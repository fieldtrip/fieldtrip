function test_bug731

% MEM 1500mb
% WALLTIME 00:10:00

% this function checks problems encountered with the neuromag digital
% trigger reading

cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/bug731/test_bug731.fif');
event       = ft_read_event(cfg.dataset);

% issue confirmed -> updated bug but will wait for Alex' input before
% proceeding
