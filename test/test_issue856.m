function test_issue856

% MEM 3gb
% WALLTIME 00:10:00
% DEPENDENCY ft_determine_coordsys ft_convert_coordsys ft_determine_units ft_convert_units ft_transform_geometry

%%

mri0 = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.mri'));
grad0 = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds'), 'senstype', 'meg');
headmodel0 = ft_read_headmodel(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.hdm'));

ft_hastoolbox('spm12', 1);

%%
% these have mixed units

mri1 = ft_determine_units(mri0);
grad1 = ft_determine_units(grad0);
headmodel1 = ft_determine_units(headmodel0);

%%
% convert to meter

target = 'm';
mri2 = ft_convert_units(mri1, target);
grad2 = ft_convert_units(grad1, target);
headmodel2 = ft_convert_units(headmodel1, target);

assert(~isequal(mri1.transform, mri2.transform));
assert(~isequal(grad1.coilpos, grad2.coilpos));
assert(~isequal(headmodel1.r, headmodel2.r));

%%
% convert to acpc

method = 0; % only approximate
mri3 = ft_convert_coordsys(mri2, 'acpc', method);
grad3 = ft_convert_coordsys(grad2, 'acpc', method);
headmodel3 = ft_convert_coordsys(headmodel2, 'acpc', method);

%%
% convert to millimeter

target = 'mm';
mri2 = ft_convert_units(mri1, target);
grad2 = ft_convert_units(grad1, target);
headmodel2 = ft_convert_units(headmodel1, target);

%%
% spm_affreg
mri4 = ft_convert_coordsys(mri2, 'acpc', 1);

% this round-trip is not exactly identical
mri4b = ft_convert_coordsys(mri4, 'ctf', 0);
assert(~isequal(mri2.transform, mri4b.transform));

%%
% smp_normalise
mri5 = ft_convert_coordsys(mri2, 'acpc', 2);

% this round-trip is not exactly identical
mri5b = ft_convert_coordsys(mri5, 'ctf', 0);
assert(~isequal(mri2.transform, mri5b.transform));
