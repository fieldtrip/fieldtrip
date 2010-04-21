function [v] = matlabversion;

% MATLABVERSION returns the Matlab version as a number
%
% Use as
%  v = matlabversion;
%
% An example for using this function is
%
%  if matlabversion<7
%    % do something specific to support an old matlab version
%  else
%    % do the default operation
%  end
%
% See also VERSION, VER

% Copyright (C) 2006, Robert Oostenveld
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

s = ver('matlab');
v = s.Version;

if ischar(v) 
  % try converting to a number
  n = str2num(v);
  if isempty(n)
    switch v
    case '6.5.1'
      n = 6.5; % this is accurate enough
    case '7.0.1'
      n = 7.0; % this is accurate enough
    case '7.0.4'
      n = 7.0; % this is accurate enough
    otherwise  
      warning('cannot convert matlab version into a number');
      v = v;
    end
  end
  v = n;
end

