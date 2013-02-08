function gradnew = ft_convert_grad(gradold, fieldstrength, distance, scaling)

% FT_CONVERT_GRAD updates the gradiometer definition to the desired
% units of distance and fieldstrength. Furthermore, it can convert the
% gradiometer channels from 'fieldstrength' to 'fieldstrength/distance', e.g.
% from T to T/m.
%
% Use as
%   grad = ft_convert_grad(grad, fieldstrength, geometry, gradiometer)
% where the input gradiometer array is according to FT_DATATYPE_SENS and
%   fieldstrength = string, can be 'T' or 'fT'
%   distance      = string, can be 'm', 'cm' or 'mm'
%   distance      = string, can be 'm', 'cm' or 'mm'
%   scaling       = string, can be 'fieldstrength' or 'fieldstrength/distance'
%
% See also FT_CHANTYPE

% Copyright (C) 2013, Robert Oostenveld
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

% perform some sanity checks
if ~ft_senstype(gradold, 'meg')
  error('unsupported type of gradiometer array "%s"', ft_senstype(gradold));
end

if ~any(strcmp(fieldstrength, {'T' 'mT' 'uT' 'nT' 'pT' 'fT'}))
  error('unsupported unit of fieldstrength "%s"', fieldstrength);
end

if ~any(strcmp(distance, {'m' 'dm' 'cm' 'mm'}))
  error('unsupported unit of distance "%s"', distance);
end

sel_m  = ~cellfun(@isempty, regexp(gradold.chanunit, '/m$'));
sel_dm = ~cellfun(@isempty, regexp(gradold.chanunit, '/dm$'));
sel_cm = ~cellfun(@isempty, regexp(gradold.chanunit, '/cm$'));
sel_mm = ~cellfun(@isempty, regexp(gradold.chanunit, '/mm$'));

if     strcmp(gradold.unit, 'm') && (any(sel_dm) || any(sel_cm) || any(sel_mm))
  error('inconsistent units in input gradiometer');
elseif strcmp(gradold.unit, 'dm') && (any(sel_m) || any(sel_cm) || any(sel_mm))
  error('inconsistent units in input gradiometer');
elseif strcmp(gradold.unit, 'cm') && (any(sel_m) || any(sel_dm) || any(sel_mm))
  error('inconsistent units in input gradiometer');
elseif strcmp(gradold.unit, 'mm') && (any(sel_m) || any(sel_dm) || any(sel_cm))
  error('inconsistent units in input gradiometer');
end

% update the units of distance
gradnew = ft_convert_units(gradold, distance);

% update the units of fieldstrength
nchan = length(gradnew.label);
for i=1:nchan
  if ~isempty(regexp(gradnew.chanunit{i}, ['/' distance '$'], 'once'))
    % this channel is expressed as fieldstrength per distance
    gradnew.tra(i,:)    = gradnew.tra(i,:) * scalingfactor(gradnew.chanunit{i}, [fieldstrength '/' distance]);
    gradnew.chanunit{i} = [fieldstrength '/' distance];
  elseif ~isempty(regexp(gradnew.chanunit{i}, 'T$', 'once'))
    % this channel is expressed as fieldstrength
    gradnew.tra(i,:)    = gradnew.tra(i,:) * scalingfactor(gradnew.chanunit{i}, fieldstrength);
    gradnew.chanunit{i} = fieldstrength;
  else
    error('unexpected channel unit "%s" in channel %d', i, gradnew.chanunit{i});
  end
end

% update the gradiometer scaling
if strcmp(scaling, 'fieldstrength')
  for i=1:nchan
    if strcmp(gradnew.chanunit{i}, [fieldstrength '/' distance])
      % this channel is expressed as fieldstrength per distance
      coil = find(abs(gradnew.tra(i,:))~=0);
      if length(coil)~=2
        error('unexpected number of coils contributing to channel %d', i);
      end
      baseline = norm(gradnew.coilpos(coil(1),:) - gradnew.coilpos(coil(2),:));
      gradnew.tra(i,:)    = gradnew.tra(i,:)*baseline;  % scale with the baseline distance
      gradnew.chanunit{i} = fieldstrength;
    elseif strcmp(gradnew.chanunit{i}, fieldstrength)
      % no conversion needed
    else
      error('unexpected channel unit "%s" in channel %d', i, gradnew.chanunit{i});
    end % if
  end % for
  
elseif strcmp(scaling, 'fieldstrength/distance')
  for i=1:nchan
    if strcmp(gradnew.chanunit{i}, fieldstrength)
      % this channel is expressed as fieldstrength
      coil = find(abs(gradnew.tra(i,:))~=0);
      if length(coil)==1
        % this is a magnetometer channel, no conversion needed
        continue
      elseif length(coil)~=2
        error('unexpected number of coils (%d) contributing to channel %s (%d)', length(coil), gradnew.label{i}, i);
      end
      baseline = norm(gradnew.coilpos(coil(1),:) - gradnew.coilpos(coil(2),:));
      gradnew.tra(i,:)    = gradnew.tra(i,:)/baseline; % scale with the baseline distance
      gradnew.chanunit{i} = [fieldstrength '/' distance];
    elseif strcmp(gradnew.chanunit{i}, [fieldstrength '/' distance])
      % no conversion needed
    else
      error('unexpected channel unit "%s" in channel %d', i, gradnew.chanunit{i});
    end % if
  end % for
  
else
  error('incorrect specification of gradiometer scaling');
end
