function test_issue2151

% WALLTIME 00:20:00
% MEM 2gb
% DEPENDENCY ft_read_atlas

% this is a test function to check whether ft_read_atlas manages to read in
% a SimNIBS based segmented volume


%path =  '/home/lau/simnibs_examples/v4/ernie/m2m_ernie';
%atlas_filenamelabels = 'final_tissues_LUT.txt';

path = '/home/common/matlab/fieldtrip/data/test/pull2156';
atlas_filename = 'final_tissues.nii.gz';

atlas = ft_read_atlas(fullfile(path, atlas_filename));
