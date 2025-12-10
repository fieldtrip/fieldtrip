function test_bug1053

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_datatype_sens
% DATA private

% the following was enough to reproduce the bug
cd(dccnpath('/project/3031000.02/test'))
load bug1053.mat

ft_datatype_sens(elec);

