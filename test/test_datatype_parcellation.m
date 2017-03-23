function test_datatype_parcellation

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_datatype_parcellation ft_datatype

% See also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1652 which
% includes an elaborate discussion to the 2012 version of the
% segmentation and parcellation structures.

clear all

% construct the data structures that are given as examples in the function help
[pnt, tri] = icosahedron162;

example1 = [];
example1.pos = pnt;
example1.tri = tri;
example1.coordsys = 'ctf';
example1.unit = 'mm';
example1.brodmann = zeros(162,1);
example1.brodmannlabel = {'Brodmann Area 1', 'Brodmann Area 2', 'Brodmann Area 3'};
example1.brodmann(1) = 1;
example1.brodmann(2) = 2;
example1.brodmann(3) = 3;

example2 = [];
example2.pos = pnt;
example2.tri = tri;
example2.coordsys = 'ctf';
example2.unit = 'mm';
example2.Brodmann_Area_1 = false(162,1);
example2.Brodmann_Area_2 = false(162,1);
example2.Brodmann_Area_3 = false(162,1);
example2.Brodmann_Area_1(1) = true;
example2.Brodmann_Area_2(2) = true;
example2.Brodmann_Area_3(3) = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform some checks on the example data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(ft_datatype(example1, 'parcellation'));
assert(ft_datatype(example2, 'parcellation'));

example1b = ft_checkdata(example1);
example2b = ft_checkdata(example2);

example1c = ft_checkdata(example1, 'parcellationstyle', 'indexed');
example2c = ft_checkdata(example2, 'parcellationstyle', 'indexed');

example1d = ft_checkdata(example1, 'parcellationstyle', 'probabilistic');
example2d = ft_checkdata(example2, 'parcellationstyle', 'probabilistic');
