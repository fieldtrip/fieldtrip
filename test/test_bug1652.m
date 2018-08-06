function test_bug1652

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_datatype_segmentation

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1652'));

load seg1.mat
% this contains TPMs
seg1_i  = ft_datatype_segmentation(seg1, 'segmentationstyle', 'indexed');
seg1_p  = ft_datatype_segmentation(seg1, 'segmentationstyle', 'probabilistic');
% seg1_ib = ft_datatype_segmentation(seg1, 'hasbrain', 'yes', 'segmentationstyle', 'indexed')  % this one is not supported due to the combination of the TPMs with the brain
seg1_pb = ft_datatype_segmentation(seg1, 'hasbrain', 'yes', 'segmentationstyle', 'probabilistic');
clear seg1*

load seg2.mat
% this one has three BEM tissues
% note that the skull and scalp are not consistent (one is cumulative, the other exclusive)
seg2_i  = ft_datatype_segmentation(seg2, 'segmentationstyle', 'indexed');
seg2_p  = ft_datatype_segmentation(seg2, 'segmentationstyle', 'probabilistic');
seg2_ib = ft_datatype_segmentation(seg2, 'hasbrain', 'yes', 'segmentationstyle', 'indexed');
seg2_pb = ft_datatype_segmentation(seg2, 'hasbrain', 'yes', 'segmentationstyle', 'probabilistic');
clear seg2*

load seg3.mat
% this contains three cumulative or overlapping BEM tissues
seg3_i  = ft_datatype_segmentation(seg3, 'segmentationstyle', 'indexed');
seg3_p  = ft_datatype_segmentation(seg3, 'segmentationstyle', 'probabilistic');
seg3_ib = ft_datatype_segmentation(seg3, 'hasbrain', 'yes', 'segmentationstyle', 'indexed');
seg3_pb = ft_datatype_segmentation(seg3, 'hasbrain', 'yes', 'segmentationstyle', 'probabilistic');
clear seg3*

load seg4.mat
% this contains three exclusive BEM tissues
seg4_i  = ft_datatype_segmentation(seg4, 'segmentationstyle', 'indexed');
seg4_p  = ft_datatype_segmentation(seg4, 'segmentationstyle', 'probabilistic');
seg4_ib = ft_datatype_segmentation(seg4, 'hasbrain', 'yes', 'segmentationstyle', 'indexed');
seg4_pb = ft_datatype_segmentation(seg4, 'hasbrain', 'yes', 'segmentationstyle', 'probabilistic');
clear seg4*

load seg5.mat
% this one is an indexed representation of the BEM tissues, without labels
seg5_i  = ft_datatype_segmentation(seg5, 'segmentationstyle', 'indexed');
seg5_p  = ft_datatype_segmentation(seg5, 'segmentationstyle', 'probabilistic');

seg5.seglabel = {'scalp', 'skull', 'brain'};
seg5_ib = ft_datatype_segmentation(seg5, 'hasbrain', 'yes', 'segmentationstyle', 'indexed')         % this fails without the labels;
seg5_pb = ft_datatype_segmentation(seg5, 'hasbrain', 'yes', 'segmentationstyle', 'probabilistic')   % this fails without the labels;
clear seg5*

load seg6.mat
% this only contains a single boolean brain
seg6_i  = ft_datatype_segmentation(seg6, 'segmentationstyle', 'indexed');
seg6_p  = ft_datatype_segmentation(seg6, 'segmentationstyle', 'probabilistic');
seg6_ib = ft_datatype_segmentation(seg6, 'hasbrain', 'yes', 'segmentationstyle', 'indexed');
seg6_pb = ft_datatype_segmentation(seg6, 'hasbrain', 'yes', 'segmentationstyle', 'probabilistic');
clear seg6*

% try it with the AFNI atlas
atlas = ft_read_atlas('TTatlas+tlrc.BRIK');
atlas_i  = ft_datatype_segmentation(atlas, 'segmentationstyle', 'indexed');
atlas_p  = ft_datatype_segmentation(atlas, 'segmentationstyle', 'probabilistic');
