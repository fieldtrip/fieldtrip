function [object] = ft_convert_units(object, target, varargin)

% FT_CONVERT_UNITS changes the geometrical dimension to the specified SI unit.
% The units of the input object is determined from the structure field
% object.unit, or is estimated based on the spatial extend of the structure,
% e.g. a volume conduction model of the head should be approximately 20 cm large.
%
% Use as
%   [output] = ft_convert_units(input, target)
%
% The following input data structures are supported
%   electrode or gradiometer array, see FT_DATATYPE_SENS
%   volume conductor, see FT_DATATYPE_HEADMODEL
%   anatomical mri, see FT_DATATYPE_VOLUME
%   segmented mri, see FT_DATATYPE_SEGMENTATION
%   source model, see FT_DATATYPE_SOURCE and FT_PREPARE_SOURCEMODEL
%
% The possible target units are 'm', 'cm ' or 'mm'. If no target units are specified,
% this function will only determine the geometrical units of the input object.
%
% See also FT_DETERMINE_UNITS, FT_DETERMINE_COORDSYS, FT_CONVERT_COORDSYS, FT_PLOT_AXES, FT_PLOT_XXX

% Copyright (C) 2005-2020, Robert Oostenveld
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

% This function consists of three parts:
%   1) determine the input units
%   2) determine the requested scaling factor to obtain the output units
%   3) try to apply the scaling to the known geometrical elements in the input object

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input options
feedback = ft_getopt(varargin, 'feedback', false);

if isstruct(object) && numel(object)>1
  % deal with a structure array
  for i=1:numel(object)
    tmp(i) = ft_convert_units(object(i), target, varargin{:});
  end
  object = tmp;
  return
elseif iscell(object) && numel(object)>1
  % deal with a cell-array
  % this might represent combined EEG, ECoG and/or MEG
  for i=1:numel(object)
    object{i} = ft_convert_units(object{i}, target, varargin{:});
  end
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine the units of the input object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
object = ft_determine_units(object);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the scaling factor from the input units to the desired ones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(object) || isfield(object, 'unit') && isequal(object.unit, target)
  % there is nothing to do
  return
end

if istrue(feedback)
  % give some information about the conversion
  fprintf('converting units from ''%s'' to ''%s''\n', object.unit, target)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the scaling factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = ft_scalingfactor(object.unit, target);

% volume conductor model
if isfield(object, 'r'), object.r = scale * object.r; end
if isfield(object, 'o'), object.o = scale * object.o; end
if isfield(object, 'bnd') && isfield(object.bnd, 'pnt')
  for i=1:length(object.bnd)
    object.bnd(i).pnt = scale * object.bnd(i).pnt;
  end
end
if isfield(object, 'bnd') && isfield(object.bnd, 'pos')
  for i=1:length(object.bnd)
    object.bnd(i).pos = scale * object.bnd(i).pos;
  end
end

% old-fashioned gradiometer array
if isfield(object, 'pnt1'), object.pnt1 = scale * object.pnt1; end
if isfield(object, 'pnt2'), object.pnt2 = scale * object.pnt2; end
if isfield(object, 'prj'),  object.prj  = scale * object.prj;  end

% gradiometer array, electrode array, head shape or dipole grid
if isfield(object, 'pnt'),        object.pnt        = scale * object.pnt;        end
if isfield(object, 'pos'),        object.pos        = scale * object.pos;        end
if isfield(object, 'chanpos'),    object.chanpos    = scale * object.chanpos;    end
if isfield(object, 'chanposorg'), object.chanposold = scale * object.chanposorg; end % pre-2016 version
if isfield(object, 'chanposold'), object.chanposold = scale * object.chanposold; end % 2016 version and later
if isfield(object, 'coilpos'),    object.coilpos    = scale * object.coilpos;    end
if isfield(object, 'elecpos'),    object.elecpos    = scale * object.elecpos;    end

% gradiometer array that combines multiple coils in one channel
if isfield(object, 'tra') && isfield(object, 'chanunit')
  % find the gradiometer channels that are expressed as unit of field strength divided by unit of distance, e.g. T/cm
  for i=1:length(object.chanunit)
    tok = tokenize(object.chanunit{i}, '/');
    if ~isempty(regexp(object.chanunit{i}, 'm$', 'once'))
      % assume that it is T/m or so
      object.tra(i,:)    = object.tra(i,:) / scale;
      object.chanunit{i} = [tok{1} '/' target];
    elseif ~isempty(regexp(object.chanunit{i}, '[T|V]$', 'once'))
      % assume that it is T or V, don't do anything
    elseif strcmp(object.chanunit{i}, 'unknown')
      % assume that it is T or V, don't do anything
    elseif strcmp(object.chanunit{i}, 'snr')
      %
    else
      ft_error('unexpected units %s', object.chanunit{i});
    end
  end % for
end % if

% fiducials
if isfield(object, 'fid') && isfield(object.fid, 'pnt'), object.fid.pnt = scale * object.fid.pnt; end
if isfield(object, 'fid') && isfield(object.fid, 'pos'), object.fid.pos = scale * object.fid.pos; end

% dipole grid
if isfield(object, 'resolution'), object.resolution = scale * object.resolution; end

% x,y,zgrid can also be 'auto'
if isfield(object, 'xgrid') && ~ischar(object.xgrid), object.xgrid = scale * object.xgrid; end
if isfield(object, 'ygrid') && ~ischar(object.ygrid), object.ygrid = scale * object.ygrid; end
if isfield(object, 'zgrid') && ~ischar(object.zgrid), object.zgrid = scale * object.zgrid; end

% anatomical MRI or functional volume
if isfield(object, 'transform')
  H = diag([scale scale scale 1]);
  object.transform = H * object.transform;
end

if isfield(object, 'transformorig')
  H = diag([scale scale scale 1]);
  object.transformorig = H * object.transformorig;
end

% remove initial and params structure if they exist
if isfield(object, 'initial') && ~strcmp(target, 'mm')
  object = rmfield(object, 'initial');
  ft_warning('Removing field "initial" because potential transformations of normalised volumes only work if geometrical values are expressed in "mm"');
end

if isfield(object, 'params') && ~strcmp(target, 'mm')
  ft_warning('Removing field "params" because potential transformations of normalised volumes only work if geometrical values are expressed in "mm"');
  object = rmfield(object, 'params');
end

% sourcemodel obtained through mne also has a orig-field with the high number of vertices
if isfield(object, 'orig')
  if isfield(object.orig, 'pnt')
    object.orig.pnt = scale * object.orig.pnt;
  end
  if isfield(object.orig, 'pos')
    object.orig.pos = scale * object.orig.pos;
  end
end

% remember the unit
object.unit = target;
