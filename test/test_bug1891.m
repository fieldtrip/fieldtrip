function test_bug1891

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivity_wpli.m
% DATA no

% use case: cross-spectra with 200 repetitions for 20*20 channel combinations, each one frequency and time bin long
data = rand(200,400);
ft_connectivity_wpli(data, 'dojack', false);
