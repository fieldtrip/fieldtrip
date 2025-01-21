function test_bug731

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY
% DATA private

% this function checks problems encountered with the neuromag digital trigger reading

cfg.dataset = dccnpath('/project/3031000.02/test/bug731/test_bug731.fif');
event       = ft_read_event(cfg.dataset);

% issue confirmed -> updated bug but will wait for Alex' input before proceeding
