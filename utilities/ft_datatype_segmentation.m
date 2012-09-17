function segmentation = ft_datatype_segmentation(segmentation, varargin)

% FT_DATATYPE_SEGMENTATION describes the FieldTrip MATLAB structure for
% segmented volumetric data and atlases.
%
% A segmentation is a volumetric description which is usually derived from
% an anatomical MRI, which describes for each voxel the tissue type. It for
% example distinguishes between white matter, grey matter, csf, skull and
% skin. It is mainly used for masking in visualization, construction of
% volume conduction models and for construction of cortical sheets. An
% volume-based atlas is basically a very detailled segmentation with
% anatomical labels for each tissue type.
%
% For example, the AFNI TTatlas+tlrc segmented brain atlas (which can be
% created with FT_PREPARE_ATLAS) looks like this
%
%              dim: [161 191 141]        the size of the 3D volume in voxels 
%        transform: [4x4 double]         affine transformation matrix for mapping the voxel coordinates to head coordinate system
%            coord: 'tal'                the transformation matrix maps the voxels into this (head) coordinate system
%             unit: 'mm'                 the units in which the coordinate system is expressed
%           brick0: [161x191x141 uint8]  integer values from 1 to N, the value 0 means unknown
%           brick1: [161x191x141 uint8]  integer values from 1 to M, the value 0 means unknown
%      brick0label: {Nx1 cell}
%      brick1label: {Mx1 cell}
%
% An example of a whole-brain anatomical MRI that was segmented using FT_VOLUMESEGMENT looks like this
%
%         dim: [256 256 256]         the size of the 3D volume in voxels
%   transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%    coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%        unit: 'mm'                  the units in which the coordinate system is expressed
%        gray: [256x256x256 double]  probabilistic map of the gray matter
%       white: [256x256x256 double]  probabilistic map of the white matter
%         csf: [256x256x256 double]  probabilistic map of the cerebrospinal fluid
%
% An example segmentation with binary values that can be used for
% construction of a BEM volume conduction model of the head looks like this
%
%           dim: [256 256 256]         the dimensionality of the 3D volume
%     transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%      coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%          unit: 'mm'                  the units in which the coordinate system is expressed
%         brain: [256x256x256 logical] binary map representing the voxels which belong to the brain  
%         scalp: [256x256x256 logical] binary map representing the voxels which belong to the scalp
%         skull: [256x256x256 logical] binary map representing the voxels which belong to the skull
%
% The only difference to the volume data representation is that the segmentation
% structure contains the additional fields XXX. See FT_DATATYPE_VOLUME for
% further details.
%