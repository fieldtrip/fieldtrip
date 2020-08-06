function [volume] = ft_datatype_volume(volume, varargin)

% FT_DATATYPE_VOLUME describes the FieldTrip MATLAB structure for volumetric data.
%
% The volume data structure represents data on a regular volumetric
% 3-D grid, like an anatomical MRI, a functional MRI, etc. It can
% also represent a source reconstructed estimate of the activity
% measured with MEG. In this case the source reconstruction is estimated
% or interpolated on the regular 3-D dipole grid (like a box).
%
% An example volume structure is
%       anatomy: [181x217x181 double]  the numeric data, in this case anatomical information
%           dim: [181 217 181]         the dimensionality of the 3D volume
%     transform: [4x4 double]          affine transformation matrix for mapping the voxel coordinates to the head coordinate system
%          unit: 'mm'                  geometrical units of the coordinate system
%      coordsys: 'ctf'                 description of the coordinate system
%
% Required fields:
%   - transform, dim
%
% Optional fields:
%   - anatomy, prob, stat, grey, white, csf, or any other field with dimensions that are consistent with dim
%   - unit, size, coordsys
%
% Deprecated fields:
%   - dimord
%
% Obsoleted fields:
%   - none
%
% Revision history:
%
% (2014) The subfields in the avg and trial fields are now present in the
% main structure, e.g. source.avg.pow is now source.pow. Furthermore, the
% inside is always represented as logical array.
%
% (2012b) Ensure that the anatomy-field (if present) does not contain
% infinite values.
%
% (2012) A placeholder 2012 version was created that ensured the axes
% of the coordinate system to be right-handed. This actually never
% has made it to the default version. An executive decision regarding
% this has not been made as far as I (JM) am aware, and probably it's
% a more principled approach to keep the handedness free, so don't mess
% with it here. However, keep this snippet of code for reference.
%
% (2011) The dimord field was deprecated and we agreed that volume
% data should be 3-dimensional and not N-dimensional with arbitary
% dimensions. In case time-frequency recolved data has to be represented
% on a 3-d grid, the source representation should be used.
%
% (2010) The dimord field was added by some functions, but not by all
%
% (2003) The initial version was defined
%
% See also FT_DATATYPE, FT_DATATYPE_DIP, FT_DATATYPE_SOURCE

% Copyright (C) 2011-2015, Robert Oostenveld
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
version = ft_getopt(varargin, 'version', 'latest');

if strcmp(version, 'latest')
  version = '2014';
end

if isempty(volume)
  return;
end

% it should  never have contained these, but they might be present due to an unclear
% distinction between the volume and the source representation
if isfield(volume, 'xgrid'),     volume = rmfield(volume, 'xgrid');     end
if isfield(volume, 'ygrid'),     volume = rmfield(volume, 'ygrid');     end
if isfield(volume, 'zgrid'),     volume = rmfield(volume, 'zgrid');     end
if isfield(volume, 'frequency'), volume = rmfield(volume, 'frequency'); end
if isfield(volume, 'latency'),   volume = rmfield(volume, 'latency');   end

if isfield(volume, 'pos')
  if ~isfield(volume, 'dim')
    volume.dim = pos2dim(volume.pos);
  end
  assert(prod(volume.dim)==size(volume.pos,1), 'dimensions are inconsistent with number of grid positions');
  if  ~isfield(volume, 'transform')
    volume.transform = pos2transform(volume.pos, volume.dim);
  end
  volume = rmfield(volume, 'pos');
end

switch version
  case '2014'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(volume, 'coordsys')
      % ensure that it is in lower case
      volume.coordsys = lower(volume.coordsys);
    end
    
    if isfield(volume, 'unit')
      % ensure that it is in lower case
      volume.unit = lower(volume.unit);
    end
    
    if isfield(volume, 'dimord')
      volume = rmfield(volume, 'dimord');
    end

    if isfield(volume, 'anatomy')
      volume.anatomy(~isfinite(volume.anatomy)) = 0;
    end

    if isfield(volume, 'avg') && isstruct(volume.avg)
      % move the average fields to the main structure
      fn = fieldnames(volume.avg);
      for i=1:length(fn)
        volume.(fn{i}) = volume.avg.(fn{i});
      end
      volume = rmfield(volume, 'avg');
    end

    % ensure that it is always logical
    volume = fixinside(volume, 'logical');

    fn = getdatfield(volume);
    for i=1:numel(fn)
      try
        volume.(fn{i}) = reshape(volume.(fn{i}), volume.dim);
      catch
        ft_notice('could not reshape "%s" to the dimensions of the volume', fn{i});
      end
    end

  case '2012b'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(volume, 'dimord')
      volume = rmfield(volume, 'dimord');
    end

    if isfield(volume, 'anatomy')
      volume.anatomy(~isfinite(volume.anatomy)) = 0;
    end

  case '2012'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % THIS ONE DOES NOT SEEM TO HAVE EVER BEEN USED
    % HOWEVER, KEEP IT FOR DOCUMENTATION PURPOSES

    if isfield(volume, 'dimord')
      volume = rmfield(volume, 'dimord');
    end

    % ensure the axes system in the transformation matrix to be
    % right-handed
    volume = volumeflip(volume, 'right');

  case '2011'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(volume, 'dimord')
      volume = rmfield(volume, 'dimord');
    end

  case '2010'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this might have been N-dimensional and contained a dimord, but in general cannot
    % be reconstructed on the fly

  case '2003'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(volume, 'dimord')
      volume = rmfield(volume, 'dimord');
    end

  otherwise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ft_error('unsupported version "%s" for volume datatype', version);
end
