function test_datatype_segmentation

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_datatype_segmentation ft_datatype

% See also http =//bugzilla.fcdonders.nl/show_bug.cgi?id=1652 which
% includes an elaborate discussion to the 2012 version of the
% segmentation and parcellation structures.

clear all

% construct the data structures that are given as examples in the function help

% For example, the AFNI TTatlas+tlrc segmented brain atlas (which can be created with FT_PREPARE_ATLAS) looks like this
example1.dim = [161 191 141];
example1.transform = eye(4);
example1.coordsys = 'tal';
example1.unit = 'mm';
example1.brick0 = zeros(161,191,141, 'uint8');
example1.brick1 = zeros(161,191,141, 'uint8');
example1.brick0label = {'a', 'b', 'c'};
example1.brick1label = {'d', 'e'};
% note that brick1 stays empty
example1.brick0(1) = 1;
example1.brick0(2) = 2;
example1.brick0(3) = 3;

% An example of a whole-brain anatomical MRI that was segmented using FT_VOLUMESEGMENT looks like this
example2.dim = [256 256 256];
example2.transform = eye(4);
example2.coordsys = 'ctf';
example2.unit = 'mm';
example2.gray  = zeros(256,256,256);
example2.white = zeros(256,256,256);
example2.csf   = zeros(256,256,256);
% the values can be between 0 and 1 (inclusive)
example2.gray (1) = 0.7;
example2.white(2) = 0.8;
example2.csf  (3) = 0.9;

% An example segmentation with binary values that can be used for construction of a BEM volume conduction model of the head looks like this
example3.dim = [256 256 256];
example3.transform = eye(4);
example3.coordsys = 'ctf';
example3.unit = 'mm';
example3.brain = false(256,256,256);
example3.scalp = false(256,256,256);
example3.skull = false(256,256,256);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform some checks on the example data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(ft_datatype(example1, 'segmentation'));
assert(ft_datatype(example2, 'segmentation'));
assert(ft_datatype(example3, 'segmentation'));

example1b = ft_checkdata(example1);
example2b = ft_checkdata(example2);
example3b = ft_checkdata(example3);

example1c = ft_checkdata(example1, 'segmentationstyle', 'indexed');
example2c = ft_checkdata(example2, 'segmentationstyle', 'indexed');
example3c = ft_checkdata(example3, 'segmentationstyle', 'indexed');

example1d = ft_checkdata(example1, 'segmentationstyle', 'probabilistic');
example2d = ft_checkdata(example2, 'segmentationstyle', 'probabilistic');
example3d = ft_checkdata(example3, 'segmentationstyle', 'probabilistic');
