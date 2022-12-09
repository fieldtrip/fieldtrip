%% INITIALISE

clear variables
restoredefaultpath

addpath /home/lau/matlab/fieldtrip
ft_defaults

%% READ ATLAS

path =  '/home/lau/simnibs_examples/v4/ernie/m2m_ernie';
atlas_filenamelabels = 'final_tissues_LUT.txt';
atlas_filenamemesh = 'final_tissues.nii.gz';

atlas = ft_read_atlas(fullfile(path, ...
                     {atlas_filenamelabels atlas_filenamemesh}));

