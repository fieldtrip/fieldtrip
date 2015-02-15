function parcellation = ft_datatype_parcellation(parcellation, varargin)

% FT_DATATYPE_PARCELLATION describes the FieldTrip MATLAB structure for parcellated
% cortex-based data and atlases. A parcellation can either be indexed or probabilistic
% (see below).
%
% A parcellation describes the tissue types for each of the surface elements.
% Parcellations are often, but not always labeled. A parcellatoin can be used to
% estimate the activity from MEG data in a known region of interest. A surface-based
% atlas is basically a very detailled parcellation with an anatomical label for each
% vertex.
%
% An example of a surface based Brodmann parcellation looks like this
%
%              pos: [8192x3]         positions of the vertices forming the cortical sheet
%              tri: [16382x3]        triangles of the cortical sheet
%         coordsys: 'ctf'            the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'             the units in which the coordinate system is expressed
%         brodmann: [8192x1 uint8]   values from 1 to N, the value 0 means unknown
%    brodmannlabel: {Nx1 cell}
%
% An alternative representation of this parcellation is
%
%              pos: [8192x3]           positions of the vertices forming the cortical sheet
%              tri: [16382x3]          triangles of the cortical sheet
%         coordsys: 'ctf'              the (head) coordinate system in which the vertex positions are expressed
%             unit: 'mm'               the units in which the coordinate system is expressed
%  Brodmann_Area_1: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  Brodmann_Area_2: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  Brodmann_Area_3: [8192x1 logical]   binary map representing the voxels belonging to the specific area
%  ...
%
% The examples above demonstrate that a parcellation can be indexed, i.e. consisting of
% subsequent integer numbers (1, 2, ...) or probabilistic, consisting of real numbers
% ranging from 0 to 1 that represent probabilities between 0% and 100%. An extreme case
% is one where the probability is either 0 or 1, in which case the probability can be
% represented as a binary or logical array.
%
% The only difference to the source data structure is that the parcellation structure
% contains the additional fields xxx and xxxlabel. See FT_DATATYPE_SOURCE for further
% details.
%
% Required fields:
%   - pos
%
% Optional fields:
%   - tri, coordsys, unit
%
% Deprecated fields:
%   - none
%
% Obsoleted fields:
%   - none
%
% Revision history:
% (2012/latest) The initial version was defined in accordance with the representation of
% a voxel-based segmentation.
%
% See also FT_DATATYPE, FT_DATATYPE_SOURCE, FT_DATATYPE_SEGMENTATION

% Copyright (C) 2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
parcellationstyle = ft_getopt(varargin, 'parcellationstyle');  % can be indexed or probabilistic

if strcmp(version, 'latest')
  parcelversion = '2012';
  sourceversion = 'latest';
  clear version
else
  parcelversion = version;
  sourceversion = version;
  clear version
end

if isempty(parcellation)
  return;
end

switch parcelversion
  case '2012'
    
    if isfield(parcellation, 'pnt')
      parcellation.pos = parcellation.pnt;
      parcellation = rmfield(parcellation, 'pnt');
    end
    
    % convert the inside/outside fields, they should be logical rather than an index
    if isfield(parcellation, 'inside')
      parcellation = fixinside(parcellation, 'logical');
    end
    
    dim = size(parcellation.pos,1);
    
    % make a list of fields that represent a parcellation
    fn = fieldnames(parcellation);
    fn = setdiff(fn, 'inside'); % exclude the inside field from any conversions
    sel = false(size(fn));
    for i=1:numel(fn)
      sel(i) = isnumeric(parcellation.(fn{i})) && numel(parcellation.(fn{i}))==dim;
    end
    % only consider numeric fields of the correct size
    fn = fn(sel);
    
    % determine whether the style of the input fields is probabilistic or indexed
    [indexed, probabilistic] = determine_segmentationstyle(parcellation, fn, dim);
    
    % ignore the fields that do not contain a parcellation
    sel = indexed | probabilistic;
    fn            = fn(sel);
    indexed       = indexed(sel);
    probabilistic = probabilistic(sel);

    if ~any(probabilistic) && ~any(indexed)
      % rather than being described with a tissue label for each vertex
      % it can also be described with a tissue label for each surface or volme element
      for i = 1:length(fn)
        fname = fn{i};
        switch fname
          case 'tri'
            dim = size(parcellation.tri,1);
          case 'hex'
            dim = size(parcellation.hex,1);
          case 'tet'
            dim = size(parcellation.tet,1);
        end
      end
      [indexed, probabilistic] = determine_segmentationstyle(parcellation, fn, dim);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ensure that the parcellation is internally consistent
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if any(probabilistic)
      parcellation = fixsegmentation(parcellation, fn(probabilistic), 'probabilistic');
    end
    
    if any(indexed)
      parcellation = fixsegmentation(parcellation, fn(indexed), 'indexed');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert the parcellation to the desired style
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(parcellationstyle, 'indexed') && any(probabilistic)
      parcellation  = convert_segmentationstyle(parcellation, fn(probabilistic), [dim 1], 'indexed');
    elseif strcmp(parcellationstyle, 'probabilistic') && any(indexed)
      parcellation  = convert_segmentationstyle(parcellation, fn(indexed), [dim 1], 'probabilistic');
    end % converting converting to desired style
    
  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    error('unsupported version "%s" for parcellation datatype', parcelversion);
end

% the parcellation is a speciat type of volume structure, so ensure that it also fulfills the requirements for that
parcellation = ft_datatype_source(parcellation, 'version', sourceversion);

% For the pass through ft_datatype_volume it is perhaps neccessary to remove
% the fields that are specific for the parcellation and add them later again.
% At this moment ft_datatype_volume nicely passes all fields, so there is no
% special handling of the parcellation fields needed.

