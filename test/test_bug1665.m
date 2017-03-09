function test_bug1665

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_volumesegment ft_volumenormalise align_ctf2spm ft_hastoolbox
% ft_checkdata

% this function tests whether align_ctf2spm works robustly when the input
% MRI contains NaNs

ft_hastoolbox('spm8',1,0);
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1665/segmentedS2.mat'));
struct_reslice = ft_checkdata(struct_reslice, 'datatype', 'volume');

cd(dccnpath('/home/common/matlab/fieldtrip/utilities/private'));
%cd(dccnpath('/home/language/jansch/matlab/fieldtrip/utilities/private'));
output = align_ctf2spm(struct_reslice);

struct_reslice.coordsys = 'itab';
output = align_itab2spm(struct_reslice);
