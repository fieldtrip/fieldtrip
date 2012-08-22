function test_bug1665

% TEST test_bug1665
% TEST ft_volumesegment ft_volumenormalise align_ctf2spm

% this function tests whether align_ctf2spm works robustly when the input
% MRI contains NaNs

load('/home/common/matlab/fieldtrip/data/test/bug1665/segmentedS2.mat');
struct_reslice = ft_checkdata(struct_reslice, 'datatype', 'volume');

cd('/home/common/matlab/fieldtrip/utilities/private');
%cd('/home/language/jansch/matlab/fieldtrip/utilities/private');
output = align_ctf2spm(struct_reslice);

struct_reslice.coordsys = 'itab';
output = align_itab2spm(struct_reslice);
