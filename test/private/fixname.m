function str = fixname(str)

% FIXNAME changes all inappropriate characters in a sting into '_'
% such that it can be used as a filename or as a structure field name. If
% the string begins with a numeric digit, an 'x' is prepended.
%
% Use as
%   str = fixname(str)
%
%
% MATLAB 2014a introduces the matlab.lang.makeValidName and
% matlab.lang.makeUniqueStrings functions for constructing unique MATLAB identifiers,
% but this particular implementation also works with older MATLAB versions.
%
% See also DEBLANK

% Copyright (C) 2012-2014, Robert Oostenveld
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

if int8(str(1))<58 && int8(str(1))>47
  % the string begins with a digit, prepend an 'x'
  str = ['x' str];
end

% truncate the string if it's too long: MATLAB maximizes the string length to 63
% characters (and throws a warning when truncating)
if numel(str)>63
  ft_warning(sprintf('%s exceeds MATLAB''s maximum name length of 63 characters and has been truncated to %s',str,str(1:63)));
  str = str(1:63);
end
