function test_bug1405

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1405
% TEST ft_checkdata ft_senstype

load /home/common/matlab/fieldtrip/template/headmodel/standard_mri.mat
mri = ft_checkdata(mri, 'hasunit', 'yes');

