function [segmentation] = ft_datatype_segmentation(segmentation, varargin)

% FT_DATATYPE_SEGMENTATION describes the FieldTrip MATLAB structure for segmented
% voxel-based data and atlasses. A segmentation can either be indexed or probabilistic
% (see below).
%
% A segmentation is a volumetric description which is usually derived from an anatomical
% MRI, which describes for each voxel the tissue type. It for example distinguishes
% between white matter, grey matter, csf, skull and skin. It is mainly used for masking
% in visualization, construction of volume conduction models and for construction of
% cortical sheets. An volume-based atlas is basically a very detailed segmentation with
% an anatomical label for each voxel.
%
% For example, the AFNI TTatlas+tlrc segmented brain atlas (which can be created
% with FT_READ_ATLAS) looks like this
%
%              dim: [161 191 141]        the size of the 3D volume in voxels
%        transform: [4x4 double]         affine transformation matrix for mapping the voxel coordinates to head coordinate system
%         coordsys: 'tal'                the transformation matrix maps the voxels into this (head) coordinate system
%             unit: 'mm'                 the units in which the coordinate system is expressed
%           brick0: [161x191x141 uint8]  integer values from 1 to N, the value 0 means unknown
%           brick1: [161x191x141 uint8]  integer values from 1 to M, the value 0 means unknown
%      brick0label: {Nx1 cell}
%      brick1label: {Mx1 cell}
%
%
% An example segmentation with binary values that can be used for construction of a
% BEM volume conduction model of the head looks like this
%
%           dim: [256 256 256]         the dimensionality of the 3D volume
%     transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%      coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%          unit: 'mm'                  the units in which the coordinate system is expressed
%         brain: [256x256x256 logical] binary map representing the voxels which belong to the brain
%         scalp: [256x256x256 logical] binary map representing the voxels which belong to the scalp
%         skull: [256x256x256 logical] binary map representing the voxels which belong to the skull
%
% An example of a whole-brain anatomical MRI that was segmented using FT_VOLUMESEGMENT
% looks like this
%
%         dim: [256 256 256]         the size of the 3D volume in voxels
%   transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to head coordinate system
%    coordsys: 'ctf'                 the transformation matrix maps the voxels into this (head) coordinate system
%        unit: 'mm'                  the units in which the coordinate system is expressed
%        gray: [256x256x256 double]  probabilistic map of the gray matter
%       white: [256x256x256 double]  probabilistic map of the white matter
%         csf: [256x256x256 double]  probabilistic map of the cerebrospinal fluid
%
% The examples above demonstrate that a segmentation can be indexed, i.e. consisting of
% subsequent integer numbers (1, 2, ...) or probabilistic, consisting of real numbers
% ranging from 0 to 1 that represent probabilities between 0% and 100%. An extreme case
% is one where the probability is either 0 or 1, in which case the probability can be
% represented as a binary or logical array.
%
% The only difference to the volume data representation is that the segmentation
% structure contains the additional fields xxx and xxxlabel. See FT_DATATYPE_VOLUME for
% further details.
%
% Required fields:
%   - dim, transform
%
% Optional fields:
%   - coordsys, unit
%
% Deprecated fields:
%   - none
%
% Obsoleted fields:
%   - none
%
% Revision history:
% (2012/latest) The explicit distunction between the indexed and probabilistic
% representation was made. For the indexed representation the additional
% xxxlabel cell-array was introduced.
%
% (2005) The initial version was defined.
%
% See also FT_DATATYPE, FT_DATATYPE_VOLUME, FT_DATATYPE_PARCELLATION

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% get the optional input arguments, which should be specified as key-value pairs
version           = ft_getopt(varargin, 'version', 'latest');
segmentationstyle = ft_getopt(varargin, 'segmentationstyle');  % can be indexed or probabilistic
hasbrain          = ft_getopt(varargin, 'hasbrain', 'no');     % no means that it is not required, if present it won't be removed

% convert from string into boolean
hasbrain = istrue(hasbrain);

if strcmp(version, 'latest')
  segversion = '2012';
  volversion = 'latest';
  clear version
else
  segversion = version;
  volversion = version;
  clear version
end

if isempty(segmentation)
  return;
end

switch segversion
  case '2012'
    % determine whether the style of the input fields is probabilistic or indexed
    fn = fieldnames(segmentation);
    fn = setdiff(fn, 'inside'); % exclude the inside field from any conversions
    [indexed, probabilistic] = determine_segmentationstyle(segmentation, fn, segmentation.dim);

    % ignore the fields that do not contain a segmentation
    sel = indexed | probabilistic;
    fn            = fn(sel);
    indexed       = indexed(sel);
    probabilistic = probabilistic(sel);

    % convert from an exclusive to cumulative representation
    % this is only only for demonstration purposes
    % for i=1:length(sel)
    %   segmentation.(fn{sel(i)}) = volumefillholes(segmentation.(fn{sel(i)}));
    % end

    [dum, i] = intersect(fn, {'scalp', 'skull', 'brain'});
    if numel(i)==3
      % put them in the preferred order
      fn(i) = {'scalp', 'skull', 'brain'};
    end
    [dum, i] = intersect(fn, {'skin', 'skull', 'brain'}); % this is not likely
    if numel(i)==3
      % put them in the preferred order
      fn(i) = {'skin', 'skull', 'brain'};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that the segmentation is internally consistent
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if any(probabilistic)
      segmentation = fixsegmentation(segmentation, fn(probabilistic), 'probabilistic');
    end
    if any(indexed)
      segmentation = fixsegmentation(segmentation, fn(indexed), 'indexed');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert the segmentation to the desired style
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(segmentationstyle)
      % keep it as it is
    elseif strcmp(segmentationstyle, 'indexed') && any(probabilistic)
      segmentation  = convert_segmentationstyle(segmentation, fn(probabilistic), segmentation.dim, 'indexed');
      indexed(probabilistic)       = true;  % these are now indexed
      probabilistic(probabilistic) = false; % these are now indexed
    elseif strcmp(segmentationstyle, 'probabilistic') && any(indexed)
      segmentation  = convert_segmentationstyle(segmentation, fn(indexed), segmentation.dim, 'probabilistic');
      probabilistic(indexed) = true;  % these are now probabilistic
      indexed(indexed)       = false; % these are now probabilistic
    end % converting between probabilistic and indexed

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add the brain if requested
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if hasbrain
      if all(indexed)
        fn = fieldnames(segmentation);
        sel = false(size(fn));
        for i=1:numel(fn)
          sel(i) = any(strcmp(fn, [fn{i} 'label']));
        end
        fn = fn(sel);

        if numel(fn)>1
          ft_error('cannot construct a brain mask on the fly; this requires a single indexed representation');
        else
          seg      = segmentation.(fn{1});
          seglabel = segmentation.([fn{1} 'label']);
          if ~any(strcmp(seglabel, 'brain'))
            threshold = 0.5;
            smooth    = 5;
            % ensure that the segmentation contains the brain mask, if not then construct it from gray+white+csf
            if length(intersect(seglabel, {'gray' 'white' 'csf'}))~=3
              ft_error('cannot construct a brain mask on the fly; this requires gray, white and csf');
            end
            gray  = seg==find(strcmp(seglabel, 'gray'));
            white = seg==find(strcmp(seglabel, 'white'));
            csf   = seg==find(strcmp(seglabel, 'csf'));
            brain = gray + white + csf;
            clear gray white csf seg
            brain = volumesmooth(brain,    smooth,    'brain');
            brain = volumethreshold(brain, threshold, 'brain');
            % store it in the output
            segmentation.brain = brain;
          end % try to construct the brain
        end

      elseif all(probabilistic)
        if ~isfield(segmentation, 'brain')
          if ~all(isfield(segmentation, {'gray' 'white' 'csf'}))
            ft_error('cannot construct a brain mask on the fly; this requires gray, white and csf');
          end
          threshold = 0.5;
          smooth    = 5;
          % ensure that the segmentation contains the brain mask, if not then construct it from gray+white+csf tissue probability maps
          gray  = segmentation.gray;
          white = segmentation.white;
          csf   = segmentation.csf;
          brain = gray + white + csf;
          clear gray white csf
          brain = volumesmooth(brain,    smooth,    'brain');
          brain = volumethreshold(brain, threshold, 'brain');
          % store it in the output
          segmentation.brain = brain;
        end
      else
        ft_error('cannot construct a brain mask on the fly; this requires a uniquely indexed or a uniquely probabilitic representation');
      end
    end % if hasbrain

  case '2005'
    % the only difference is that the indexed representation for xxx did not have the xxxlabel field prior to the 2012 version
    fn = fieldnames(segmentation);
    sel = ~cellfun(@isempty, regexp(fn, 'label$'));
    segmentation = rmfield(segmentation, fn(sel));
    % furthermore it corresponds to the oldest version of the volume representation
    volversion = '2003';

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ft_error('unsupported version "%s" for segmentation datatype', segversion);
end

% the segmentation is a special type of volume structure, so ensure that it also fulfills the requirements for that
segmentation = ft_datatype_volume(segmentation, 'version', volversion);

% For the pass through ft_datatype_volume it is perhaps necessary to remove
% the fields that are specific for the segmentation and add them later again.
% At this moment ft_datatype_volume nicely passes all fields, so there is no
% special handling of the segmentation fields needed.

