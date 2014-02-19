function test_bug1309

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1309

% This tests for the reliability of the ft_convert_units function when dealing
% with different input types (vol,sens,etc.)

ft_defaults;

% disable verbose output
global ft_default;
ft_default.feedback = 'no';

% artificially create timelock data
timelock_data.fsamepl = 500;
timelock_data.dimord = 'chan_time';
timelock_data.time = zeros(1, 500);
timelock_data.avg = randn(31, 500);
for i = 1:31
    timelock_data.label{i} = sprintf('label%d', i);
end

cfg = [];
cfg.elec = ft_read_sens('standard_1020.elc');
cfg.vol = ft_read_vol('standard_bem.mat');
cfg.grid.resolution = 3;
cfg.model = 'regional';
cfg.latency = [0.14,0.17];
cfg.numdipoles = 1;
timelock_data.label = cfg.elec.label(1:31);

dip = ft_dipolefitting(cfg, timelock_data);
