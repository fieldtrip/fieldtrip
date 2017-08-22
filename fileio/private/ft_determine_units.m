function [data] = ft_determine_units(data)

% FT_DETERMINE_UNITS tries to determine the units of a geometrical object by
% looking at its size and by relating this to the approximate size of the
% human head according to the following table:
%   from  0.050 to   0.500 -> meter
%   from  0.500 to   5.000 -> decimeter
%   from  5.000 to  50.000 -> centimeter
%   from 50.000 to 500.000 -> millimeter
%
% Use as
%   dataout = ft_determine_units(datain)
% where the input data structure can be
%  - an anatomical MRI
%  - an electrode or gradiometer definition
%  - a volume conduction model of the head
% or most other FieldTrip structures that represent geometrical information.
%
% See also FT_CONVERT_UNITS, FT_DETERMINE_COODSYS, FT_CONVERT_COORDSYS

% Copyright (C) 2005-2017, Robert Oostenveld
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

data  = ft_checkdata(data);
dtype = ft_datatype(data);

if isfield(data, 'unit') && ~isempty(data.unit)
  % use the units specified in the object
  unit = data.unit;
  
elseif isfield(data, 'bnd') && isfield(data.bnd, 'unit')
  
  unit = unique({data.bnd.unit});
  if ~all(strcmp(unit, unit{1}))
    ft_error('inconsistent units in the individual boundaries');
  else
    unit = unit{1};
  end
  
  % keep one representation of the units rather than keeping it with each boundary
  % the units will be reassigned further down
  data.bnd = rmfield(data.bnd, 'unit');
  
else
  % try to determine the units by looking at the size of the object
  if isfield(data, 'chanpos') && ~isempty(data.chanpos)
    siz = norm(idrange(data.chanpos));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'elecpos') && ~isempty(data.elecpos)
    siz = norm(idrange(data.elecpos));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'coilpos') && ~isempty(data.coilpos)
    siz = norm(idrange(data.coilpos));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'pnt') && ~isempty(cat(1, data.pnt))
    % the data can be a struct.array, hence the concatenation
    siz = norm(idrange(cat(1, data.pnt)));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'pos') && ~isempty(cat(1, data.pos))
    % the data can be a struct.array, hence the concatenation
    siz = norm(idrange(cat(1, data.pos)));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'transform') && ~isempty(data.transform)
    % construct the corner points of the volume in voxel and in head coordinates
    [pos_voxel, pos_head] = cornerpoints(data.dim, data.transform);
    siz = norm(idrange(pos_head));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'fid') && isfield(data.fid, 'pnt') && ~isempty(data.fid.pnt)
    siz = norm(idrange(data.fid.pnt));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'fid') && isfield(data.fid, 'pos') && ~isempty(data.fid.pos)
    siz = norm(idrange(data.fid.pos));
    unit = ft_estimate_units(siz);
    
  elseif ft_voltype(data, 'infinite')
    % this is an infinite medium volume conductor, which does not care about units
    unit = 'm';
    
  elseif ft_voltype(data,'singlesphere')
    siz = data.r;
    unit = ft_estimate_units(siz);
    
  elseif ft_voltype(data,'localspheres')
    siz = median(data.r);
    unit = ft_estimate_units(siz);
    
  elseif ft_voltype(data,'concentricspheres')
    siz = max(data.r);
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'bnd') && isstruct(data.bnd) && isfield(data.bnd(1), 'pnt') && ~isempty(data.bnd(1).pnt)
    siz = norm(idrange(data.bnd(1).pnt));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'bnd') && isstruct(data.bnd) && isfield(data.bnd(1), 'pos') && ~isempty(data.bnd(1).pos)
    siz = norm(idrange(data.bnd(1).pos));
    unit = ft_estimate_units(siz);
    
  elseif isfield(data, 'nas') && isfield(data, 'lpa') && isfield(data, 'rpa')
    pnt = [data.nas; data.lpa; data.rpa];
    siz = norm(idrange(pnt));
    unit = ft_estimate_units(siz);
    
  else
    ft_error('cannot determine geometrical units');
    
  end % recognized type of volume conduction model or sensor array
end % determine input units

% add the units to the output object
for i=1:numel(data)
  % data can be a struct-array, hence the loop
  data(i).unit = unit;
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

