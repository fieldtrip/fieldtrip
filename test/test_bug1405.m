function test_bug1405

% WALLTIME 00:03:01

% TEST test_bug1405
% TEST ft_checkdata ft_senstype

load /home/common/matlab/fieldtrip/template/headmodel/standard_mri.mat
mri = ft_checkdata(mri, 'hasunits', 'yes');

