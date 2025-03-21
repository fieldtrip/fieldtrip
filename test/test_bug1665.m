function test_bug1665

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_volumesegment ft_volumenormalise align_ctf2acpc ft_hastoolbox
% DATA private
% ft_checkdata

% this function tests whether align_ctf2acpc works robustly when the input
% MRI contains NaNs

% get the FieldTrip version and path
[ftver, ftpath] = ft_version;

ft_hastoolbox('spm8',1,0);
load(dccnpath('/project/3031000.02/test/bug1665/segmentedS2.mat'));
struct_reslice = ft_checkdata(struct_reslice, 'datatype', 'volume');

cd(fullfile(ftpath, 'utilities', 'private'));
output = align_ctf2acpc(struct_reslice);

struct_reslice.coordsys = 'itab';
output = align_neuromag2acpc(struct_reslice);
