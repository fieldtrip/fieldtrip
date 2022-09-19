function test_bug1053

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_datatype_sens

% the following was enough to reproduce the bug
cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1053.mat

ft_datatype_sens(elec);

