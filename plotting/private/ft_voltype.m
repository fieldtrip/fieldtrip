function [type] = ft_voltype(vol, desired)

% FT_VOLTYPE determines the type of volume conduction model of the head
%
% Use as
%   [type] = ft_voltype(vol)
% to get a string describing the type, or
%   [flag] = ft_voltype(vol, desired)
% to get a boolean value.
%
% For EEG the following volume conduction models are recognized
%   singlesphere       analytical single sphere model
%   concentricspheres  analytical concentric sphere model with up to 4 spheres
%   halfspace          infinite homogenous medium on one side, vacuum on the other
%   openmeeg           boundary element method, based on the OpenMEEG software
%   bemcp              boundary element method, based on the implementation from Christophe Phillips
%   dipoli             boundary element method, based on the implementation from Thom Oostendorp
%   asa                boundary element method, based on the (commercial) ASA software
%   simbio             finite element method, based on the SimBio software
%   fns                finite difference method, based on the FNS software
%   interpolate        interpolate the potential based on pre-computed leadfields
%
% and for MEG the following volume conduction models are recognized
%   singlesphere       analytical single sphere model
%   localspheres       local spheres model for MEG, one sphere per channel
%   singleshell        realisically shaped single shell approximation, based on the implementation from Guido Nolte
%   infinite           magnetic dipole in an infinite vacuum
%   interpolate        interpolate the potential based on pre-computed leadfields
%
% See also FT_COMPUTE_LEADFIELD, FT_READ_VOL, FT_HEADMODEL_BEMCP,
% FT_HEADMODEL_ASA, FT_HEADMODEL_DIPOLI, FT_HEADMODEL_SIMBIO,
% FT_HEADMODEL_FNS, FT_HEADMODEL_HALFSPACE, FT_HEADMODEL_INFINITE,
% FT_HEADMODEL_OPENMEEG, FT_HEADMODEL_SINGLESPHERE,
% FT_HEADMODEL_CONCENTRICSPHERES, FT_HEADMODEL_LOCALSPHERES,
% FT_HEADMODEL_SINGLESHELL, FT_HEADMODEL_INTERPOLATE

% Copyright (C) 2007-2012, Robert Oostenveld
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
  type = 'localspheres';
  
elseif isfield(vol, 'r') && isfield(vol, 'o') && size(vol.r,1)==size(vol.o,1) && size(vol.r,1)>4
  % this is after the spheres have been assigned to the coils
  % note that this one is easy to confuse with the concentric one
  type = 'localspheres';
  
elseif isfield(vol, 'r') && numel(vol.r)>=2 && ~isfield(vol, 'label')
  type = 'concentricspheres';
  
elseif isfield(vol, 'bnd') && isfield(vol, 'mat')
  type = 'bem'; % it could be dipoli, asa, bemcp or openmeeg
  
elseif isfield(vol, 'bnd') && isfield(vol, 'forwpar')
  type = 'singleshell';
  
elseif isfield(vol, 'bnd') && numel(vol.bnd)==1
  type = 'singleshell'; 
  
elseif isempty(vol) || (isstruct(vol) && isequal(fieldnames(vol), {'unit'}))
  % it is empty, or only contains a specification of geometrical units
  type = 'infinite';
  
else
  type = 'unknown';
  
end % if isfield(vol, 'type')

if ~isempty(desired)
  % return a boolean flag
  switch desired
    case 'bem'
      type = any(strcmp(type, {'bem', 'dipoli', 'asa', 'bemcp', 'openmeeg'}));
    otherwise
      type = any(strcmp(type, desired));
  end % switch desired
end % determine the correspondence to the desired type

% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout  = {type};
previous_argin  = current_argin;
previous_argout = current_argout;

return % voltype main()
