function varargout = ft_debug(varargin)

% FT_DEBUG prints a debug message on screen, depending on the verbosity 
% settings of the calling high-level FieldTrip function.
%
% Use as
%   ft_debug(...)
% with arguments similar to fprintf, or
%   ft_debug(msgId, ...)
% with arguments similar to warning.
%
% You can switch of all messages using
%   ft_debug off
% or for specific ones using
%   ft_debug off msgId
%
% To switch them back on, you would use 
%   ft_debug on
% or for specific ones using
%   ft_debug on msgId
% 
% Messages are only printed once per timeout period using
%   ft_debug timeout 60
%   ft_debug once
% or for specific ones using
%   ft_debug once msgId
%
% You can see the most recent messages and identifier using
%   ft_debug last
%
% You can query the current on/off/once state for all messages using
%   ft_debug query
%
% See also FT_ERROR, FT_WARNING, FT_NOTICE, FT_INFO, FT_DEBUG, ERROR, WARNING

% Copyright (C) 2017, Robert Oostenveld
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

if nargout
  varargout{:} = ft_notification(varargin{:});
elseif isequal(varargin, {'last'})
  % return an answer anyway
  varargout{1} = ft_notification(varargin{:});
else
  ft_notification(varargin{:});
end

