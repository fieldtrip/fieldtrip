function [type] = ft_voltype(vol, desired)

% FT_VOLTYPE determines the type of volume conduction model
%
% Use as
%   [type] = ft_voltype(vol)
% to get a string describing the type, or
%   [flag] = ft_voltype(vol, desired)
% to get a boolean value.
%
% See also FT_READ_VOL, FT_COMPUTE_LEADFIELD

% Copyright (C) 2007-2008, Robert Oostenveld
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
% $Id: ft_voltype.m 946 2010-04-21 17:51:16Z roboos $

% these are for remembering the type on subsequent calls with the same input arguments
persistent previous_argin previous_argout

if iscell(vol) && numel(vol)<4
  % this might represent combined EEG, ECoG and/or MEG
  type = cell(size(vol));
  if nargin<2
    desired = cell(size(vol)); % empty elements
  end
  for i=1:numel(vol)
    type{i} = ft_voltype(vol{i}, desired{i});
  end
  return
end

if nargin<2
  % ensure that all input arguments are defined
  desired = [];
end

current_argin = {vol, desired};
if isequal(current_argin, previous_argin)
  % don't do the type detection again, but return the previous values from cache
  type = previous_argout{1};
  return
end

if isfield(vol, 'type')
  % preferably the structure specifies its own type
  type = vol.type;

elseif isfield(vol, 'r') && numel(vol.r)==1 && ~isfield(vol, 'label')
  type = 'singlesphere';

elseif isfield(vol, 'r') && isfield(vol, 'o') && isfield(vol, 'label')
  % this is before the spheres have been assigned to the coils
  % and every sphere is still associated with a channel
  type = 'multisphere';

elseif isfield(vol, 'r') && isfield(vol, 'o') && size(vol.r,1)==size(vol.o,1) && size(vol.r,1)>4
  % this is after the spheres have been assigned to the coils
  % note that this one is easy to confuse with the concentric one
  type = 'multisphere';

elseif isfield(vol, 'r') && numel(vol.r)>=2 && ~isfield(vol, 'label')
  type = 'concentric';

elseif isfield(vol, 'bnd')
  type = 'bem';

elseif isempty(vol)
  type = 'infinite';

else
  type = 'unknown';

end % if isfield(vol, 'type')

if ~isempty(desired)
  % return a boolean flag
  switch desired
    case 'bem'
      type = any(strcmp(type, {'bem', 'dipoli', 'asa', 'avo', 'bemcp', 'openmeeg'}));
    otherwise
      type = any(strcmp(type, desired));
  end
end

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {type};
previous_argin  = current_argin;
previous_argout = current_argout;

return % voltype main()
