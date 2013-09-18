function str = fixname(str)

% FIXNAME changes all inappropriate characters in a sting into '_'
% such that it can be used as a filename or as a structure field name. If
% the string begins with a numeric digit, an 'x' is prepended.
%
% Use as
%   str = fixname(str)
%
% See also DEBLANK

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

str = lower(str);
str(regexp(str,'\W')) = '_';

while(str(1) == '_'),   str = str(2:end); end;   % remove all underscore at the begin of the string
while(str(end) == '_'), str = str(1:end-1); end; % remove all underscore at the end of the string

if ~isempty(str2num(str(1))) && ~isequal(str(1), 'i')
  % the string begins with a digit, prepend an 'x'
  str = ['x' str];
end
