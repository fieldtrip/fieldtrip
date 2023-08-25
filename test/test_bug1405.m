function test_bug1405

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_checkdata ft_senstype
% DATA no

[ftver, ftpath] = ft_version;
templatedir  = fullfile(ftpath, 'template');

load(fullfile(templatedir, 'headmodel', 'standard_mri.mat'));
mri = ft_checkdata(mri, 'hasunit', 'yes');

