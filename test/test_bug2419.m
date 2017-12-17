function test_bug2419

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_read_headshape

% reported bug: some particular fif-mesh files don't contain triangle area
% (use_tri_area) information, causing ft_read_headshape to crash (since it
% assumes this info to be in the file)
%
% the crash can indeed be reproduced with the file supplied. the suggested 
% fix is to check for the presence of the field.
% NOTE: the mesh in the file has an empty usetri field, suggesting that the
% mesh only contains a downsampling of the vertices, without a
% re-triangulation. this may explain the absence of the area information.

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2419.fif');
try
  bnd = ft_read_headshape(filename, 'format', 'mne_source');
catch,
  error('reading in of the headshape fif-file failed');
end
