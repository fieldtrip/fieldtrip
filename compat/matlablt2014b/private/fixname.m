function str = fixname(str, version)

% FIXNAME changes all inappropriate characters in a sting into '_'
% so that it can be used as a filename or as a field name in a structure.
% If the string begins with a digit, an 'x' is prepended.
%
% Use as
%   str = fixname(str)
%
% MATLAB 2014a introduces the matlab.lang.makeValidName and
% matlab.lang.makeUniqueStrings functions for constructing unique
% identifiers, but this particular implementation also works with
% older MATLAB versions.
%
% See also DEBLANK

% Copyright (C) 2012-2017, Robert Oostenveld
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

if nargin<2
  version = 'default';
end

switch version
  case 'default'
    if isempty(str)
      str='x';
    end
    str = lower(str);
    str(regexp(str,'\W')) = '_';
    while(str(1) == '_'),   str = str(2:end); end;   % remove all underscore at the begin of the string
    while(str(end) == '_'), str = str(1:end-1); end; % remove all underscore at the end of the string
    if int8(str(1))<58 && int8(str(1))>47
      % the string begins with a digit, prepend an 'x'
      str = ['x' str];
    end
    % truncate the string if it's too long: MATLAB maximizes the string length to 63 characters (and throws a warning when truncating)
    if numel(str)>63
      str = str(1:63);
      ft_warning(sprintf('%s exceeds MATLAB''s maximum name length of 63 characters and has been truncated to %s', str, str(1:63)));
    end
    
  case '2014a'
    str = matlab.lang.makeValidName(str);
    
  case 'X_base64encode_X'
    % this uses some code from Mathworks file exchange
    ft_hastoolbox('fileexchange', 1);
    % the following is an encoding that can be reverted
    str = strtrim(base64encode(str));
    str(str=='=') = '_';  % replace the '=' sign with '_'
    str = ['X' str 'X'];  % start and end with an 'X'
    
  case 'X_base64decode_X'
    % this uses some code from Mathworks file exchange
    ft_hastoolbox('fileexchange', 1);
    % revert the encoding
    str = str(2:end-1);   % it starts and ends with 'X'
    str(str=='_') = '=';  % the '=' sign has been replaced with '_'
    str = char(base64decode(str));
    
  otherwise
    error('unsupported version "%s"', version);
end
