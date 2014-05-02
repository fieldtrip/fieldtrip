function test_bug2556


% MEM 2gb
% WALLTIME 00:10:00

% TEST test_bug2556
% TEST ft_sourceparcellate ft_checkdata

ftpath = fileparts(which('ft_defaults'));
filename = fullfile(ftpath, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii');

aal = ft_read_atlas(filename);

aal = ft_checkdata(aal, 'datatype', 'parcellation', 'parcellationstyle', 'indexed');


