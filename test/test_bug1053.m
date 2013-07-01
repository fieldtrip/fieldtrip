function test_bug1053

% TEST test_bug1053
% TEST ft_datatype_sens

% the following was enough to reproduce the bug
cd(dccnfilename('/home/common/matlab/fieldtrip/data/test'))
load bug1053.mat

ft_datatype_sens(elec);

