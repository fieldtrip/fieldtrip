function [obj] = ft_convert_coordsys(obj, target)

% FT_CONVERT_COORDSYS changes the coordinate system of the input object to
% the specified coordinate system. The coordinate system of the input
% object is determined from the structure field object.coordsys, or need to
% be determined interactively by the user.
%
% Use as
%   [object] = ft_convert_coordsys(object, target)
%
% The following input objects are supported
%   anatomical mri
%   (not yet) electrode definition
%   (not yet) gradiometer array definition
%   (not yet) volume conductor definition
%   (not yet) dipole grid definition
%
% Possible target coordinate systems are'spm'.
%
% Note that the conversion will be an automatic one which means that it
% will be an approximate conversion, not taking into account differences in
% individual anatomies/differences in conventions where to put the
% fiducials.
%
% See also FT_DETERMINE_COORDSYS

% Copyright (C) 2005-2008, Robert Oostenveld
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
% $Id: ft_convert_coordsys.m 2885 $

if ~isfield(obj, 'coordsys') || isempty(obj.coordsys)
  obj = ft_determine_coordsys(obj, 'interactive', 'yes');
end

if strcmpi(target, obj.coordsys)
  return;
else
  switch target
    case {'spm' 'mni' 'tal'}
      switch obj.coordsys
        case {'ctf' 'bti' '4d'}
          fprintf('Converting the coordinate system from %s to %s\n', obj.coordsys, target);
          obj = align_ctf2spm(obj);
        case {'itab' 'neuromag'}
          fprintf('Converting the coordinate system from %s to %s\n', obj.coordsys, target);
          obj = align_itab2spm(obj);
        otherwise
      end %switch obj.coordsys
    otherwise
      error('conversion from %s to %s is not yet supported', obj.coordsys, target);       
  end %switch target
end
