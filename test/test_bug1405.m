function test_bug1405

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_checkdata ft_senstype

load(dccnpath('/home/common/matlab/fieldtrip/template/headmodel/standard_mri.mat'));
mri = ft_checkdata(mri, 'hasunit', 'yes');

