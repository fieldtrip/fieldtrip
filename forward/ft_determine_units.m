function [obj] = ft_determine_units(obj)

% FT_DETERMINE_UNITS tries to determine the units of a geometrical object by
% looking at its size and by relating this to the approximate size of the
% human head according to the following table:
%   from  0.050 to   0.500 -> meter
%   from  0.500 to   5.000 -> decimeter
%   from  5.000 to  50.000 -> centimeter
%   from 50.000 to 500.000 -> millimeter
%
% Use as
%   [output] = ft_determine_units(input)
%
% The following input data structures are supported
%   electrode or gradiometer array, see FT_DATATYPE_SENS
%   volume conduction model, see FT_DATATYPE_HEADMODEL
%   source model, see FT_DATATYPE_SOURCE and FT_PREPARE_SOURCEMODEL
%   anatomical mri, see FT_DATATYPE_VOLUME
%   segmented mri, see FT_DATATYPE_SEGMENTATION
%   anatomical or functional atlas, see FT_READ_ATLAS
%
% This function will add the field 'unit' to the output data structure with the
% possible values 'm', 'cm ' or 'mm'.
%
% See also FT_CONVERT_UNITS, FT_DETERMINE_COODSYS, FT_CONVERT_COORDSYS, FT_PLOT_AXES, FT_PLOT_XXX

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(obj) && numel(obj)>1
  % deal with a structure array
  for i=1:numel(obj)
    tmp(i) = ft_determine_units(obj(i));
  end
  obj = tmp;
  return
elseif iscell(obj) && numel(obj)>1
  % deal with a cell-array
  % this might represent combined EEG, ECoG and/or MEG
  for i=1:numel(obj)
    obj{i} = ft_determine_units(obj{i});
  end
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(obj, 'unit') && ~isempty(obj.unit)
  % use the units specified in the object
  unit = obj.unit;
  
elseif isfield(obj, 'bnd') && isfield(obj.bnd, 'unit')
  
  unit = unique({obj.bnd.unit});
  if ~all(strcmp(unit, unit{1}))
    ft_error('inconsistent units in the individual boundaries');
  else
    unit = unit{1};
  end
  
  % keep one representation of the units rather than keeping it with each boundary
  % the units will be reassigned further down
  obj.bnd = rmfield(obj.bnd, 'unit');
  
else
  % try to determine the units by looking at the size of the object
  if isfield(obj, 'chanpos') && ~isempty(obj.chanpos) && ~all(isnan(obj.chanpos(:)))
    siz = norm(idrange(obj.chanpos));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'elecpos') && ~isempty(obj.elecpos) && ~all(isnan(obj.elecpos(:)))
    siz = norm(idrange(obj.elecpos));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'coilpos') && ~isempty(obj.coilpos) && ~all(isnan(obj.coilpos(:)))
    siz = norm(idrange(obj.coilpos));
    unit = ft_estimate_units(siz);

  elseif isfield(obj, 'optopos') && ~isempty(obj.optopos) && ~all(isnan(obj.optopos(:)))
    siz = norm(idrange(obj.optopos));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'pnt') && ~isempty(cat(1, obj.pnt))
    % the obj can be a struct.array, hence the concatenation
    siz = norm(idrange(cat(1, obj.pnt)));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'pos') && ~isempty(cat(1, obj.pos))
    % the obj can be a struct.array, hence the concatenation
    siz = norm(idrange(cat(1, obj.pos)));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'transform') && ~isempty(obj.transform)
    % construct the corner points of the volume in voxel and in head coordinates
    [pos_voxel, pos_head] = cornerpoints(obj.dim, obj.transform);
    siz = norm(idrange(pos_head));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'fid') && isfield(obj.fid, 'pnt') && ~isempty(obj.fid.pnt)
    siz = norm(idrange(obj.fid.pnt));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'fid') && isfield(obj.fid, 'pos') && ~isempty(obj.fid.pos)
    siz = norm(idrange(obj.fid.pos));
    unit = ft_estimate_units(siz);
    
  elseif ft_headmodeltype(obj, 'infinite')
    % this is an infinite medium volume conductor, which does not care about units
    unit = 'm';
    
  elseif ft_headmodeltype(obj,'singlesphere')
    siz = obj.r;
    unit = ft_estimate_units(siz);
    
  elseif ft_headmodeltype(obj,'localspheres')
    siz = median(obj.r);
    unit = ft_estimate_units(siz);
    
  elseif ft_headmodeltype(obj,'concentricspheres')
    siz = max(obj.r);
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'bnd') && isstruct(obj.bnd) && isfield(obj.bnd(1), 'pnt') && ~isempty(obj.bnd(1).pnt)
    siz = norm(idrange(obj.bnd(1).pnt));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'bnd') && isstruct(obj.bnd) && isfield(obj.bnd(1), 'pos') && ~isempty(obj.bnd(1).pos)
    siz = norm(idrange(obj.bnd(1).pos));
    unit = ft_estimate_units(siz);
    
  elseif isfield(obj, 'nas') && isfield(obj, 'lpa') && isfield(obj, 'rpa')
    pnt = [obj.nas; obj.lpa; obj.rpa];
    siz = norm(idrange(pnt));
    unit = ft_estimate_units(siz);
    
  else
    ft_error('cannot determine geometrical units');
    
  end % recognized type of volume conduction model or sensor array
end % determine input units

% add the units to the output object
for i=1:numel(obj)
  % obj can be a struct-array, hence the loop
  obj(i).unit = unit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDRANGE interdecile range for more robust range estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = idrange(x)
keeprow=true(size(x,1),1);
for l=1:size(x,2)
  keeprow = keeprow & isfinite(x(:,l));
end
sx = sort(x(keeprow,:), 1);
ii = round(interp1([0, 1], [1, size(x(keeprow,:), 1)], [.1, .9]));  % indices for 10 & 90 percentile
r = diff(sx(ii, :));
