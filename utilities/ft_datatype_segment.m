function atlas = ft_datatype_segmentation(atlas, varargin)

% FT_DATATYPE_SEGMENT describes the FieldTrip MATLAB structure for segmented
% and parcelated data and for volumetric and cortex-based atlases.
%
% A segmentation is a volumetric description derived from an anatomical
% MRI, which describes for each voxel the tissue type.  Usually it
% separates between white matter, grey matter, csf, skull and skin.
% The white, grey and csf can also be combined into one brain
% compartment. It is mainly used for brain masking, construction of
% volume conduction models and construction of cortical sheets.
%
% A parcelation also describes the tissue types, but with much
% greater level of detail. It is for example used to represent a a
% volumetric or surface based atlas, which optionally can be determined
% for the individual subject. Parcelations are often, but not always
% labeled.
%
% For example, the AFNI TTatlas+tlrc parcelated brain atlas (which can be
% created with FT_PREPARE_ATLAS) looks like this
%
% atlas =
%              dim: [161 191 141]
%        transform: [4x4 double]
%            coord: 'tal'                the transformation matrix maps the voxels into this (head) coordinate system
%             unit: 'mm'                 the units in which the coordinate system is expressed
%           brick0: [161x191x141 uint8]  values from 1 to N, the value 0 means unknown
%           brick1: [161x191x141 uint8]  values from 1 to M, the value 0 means unknown
%      brick0label: {Nx1 cell}
%      brick1label: {Mx1 cell}
%
% An example of a surface based Brodmann parcelation looks like this
%
%              pos: [8192x3]         vertices of the cortical sheet
%              tri: [16382x3]        triangles of te cortical sheet
%            coord: 'ctf'            the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'             the units in which the coordinate system is expressed
%         brodmann: [8192x1 uint8]   values from 1 to N, the value 0 means unknown
%    brodmannlabel: {Nx1 cell}
%
% An example of a whole-brain anatomical MRI that was segmented using FT_VOLUMESEGMENT looks like this
%
%         dim: [256 256 256]         the dimensionality of the 3D volume  
%   transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%    coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%        unit: 'mm'                  the units in which the coordinate system is expressed
%        gray: [256x256x256 double]  probabilistic or binary map of the gray matter
%       white: [256x256x256 double]  probabilistic or binary map of the white matter
%         csf: [256x256x256 double]  probabilistic or binary map of the cerebrospinal fluid
%
% Required fields:
%   - to be discussed
%
% Optional fields:
%   - dim, transform, unit, coordsys
%
% Deprecated fields:
%   - descr, brick0, brick1
%
% Obsoleted fields:
%   - none
%
% Revision history:
%
% (2012/latest -> THIS IS ONLY A PROPOSAL) The representation of segmentations, parcelations and atlases has been merged.
%
% (2005) The first implementation of an atlas representation was made, based on the AFNI brick0 and brick1 representation.
%
% See also FT_DATATYPE_SOURCE, FT_DATAYPE_VOLUME

% get the optional input arguments, which should be specified as key-value pairs
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2012';
end

if isempty(atlas)
  return;
end

switch(version)
  case '2012'
    if all(isfield(atlas, {'dim', 'brick0', 'brick1', 'descr'}))
      % convert from 2005 to 2012 atlas
      
      sel0   = (atlas.descr.brick==0);
      label0 = atlas.descr.name(sel0);
      value0 = atlas.descr.value(sel0);
      % construct a new array with parcel or atlas values
      if numel(label0)<=intmax('uint8')
        new_brick0 = zeros(atlas.dim, 'uint8');
      elseif numel(label0)<=intmax('uint16')
        new_brick0 = zeros(atlas.dim, 'uint16');
      elseif numel(label0)<=intmax('uint32')
        new_brick0 = zeros(atlas.dim, 'uint32');
      else
        new_brick0 = zeros(atlas.dim);
      end
      for i=1:numel(label0)
        % replace the original value with numbers from 1 to N
        new_brick0(atlas.brick0==value0(i)) = i;
      end
      
      sel1   = (atlas.descr.brick==1);
      label1 = atlas.descr.name(sel1);
      value1 = atlas.descr.value(sel1);
      % construct a new array with parcel or atlas values
      if numel(label1)<=intmax('uint8')
        new_brick1 = zeros(atlas.dim, 'uint8');
      elseif numel(label1)<=intmax('uint16')
        new_brick1 = zeros(atlas.dim, 'uint16');
      elseif numel(label1)<=intmax('uint32')
        new_brick1 = zeros(atlas.dim, 'uint32');
      else
        new_brick1 = zeros(atlas.dim);
      end
      for i=1:numel(label1)
        % replace the original value with numbers from 1 to N
        new_brick1(atlas.brick1==value1(i)) = i;
      end
      
      atlas = rmfield(atlas, {'brick0', 'brick1', 'descr'});
      if numel(label0)>0
        % it being zero should never happen
        atlas.brick0 = new_brick0;
        atlas.brick0_label = label0;
      end
      if numel(label1)>0
        % it being zero can happen for WFU atlasses
        atlas.brick1 = new_brick1;
        atlas.brick1_label = label1;
      end
    end % conversion from 2005 to 2012
    
  case '2005'
    % these look like this, also if it contains a WFU atlas instead of the AFNI atlas.
    %
    %           dim: [161 191 141]
    %           hdr: [1x1 struct]
    %     transform: [4x4 double]
    %        brick0: [161x191x141 double]
    %        brick1: [161x191x141 double]
    %         coord: 'tal'
    %         descr: [1x1 struct]
    %           cfg: [1x1 struct]
    error('converting to the old 2005 repreaentation is not supported')
    
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for source datatype', version);
end

