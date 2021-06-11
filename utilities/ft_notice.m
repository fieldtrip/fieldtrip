function varargout = ft_notice(varargin)

% FT_NOTICE prints a notice message on screen, depending on the verbosity 
% settings of the calling high-level FieldTrip function.
%
% Use as
%   ft_notice(...)
% with arguments similar to fprintf, or
%   ft_notice(msgId, ...)
% with arguments similar to warning.
%
% You can switch of all messages using
%   ft_notice off
% or for specific ones using
%   ft_notice off msgId
%
% To switch them back on, you would use 
%   ft_notice on
% or for specific ones using
%   ft_notice on msgId
% 
% Messages are only printed once per timeout period using
%   ft_notice timeout 60
%   ft_notice once
% or for specific ones using
%   ft_notice once msgId
%
% You can see the most recent messages and identifier using
%   ft_notice last
%
% You can query the current on/off/once state for all messages using
%   ft_notice query
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
  [varargout{1:nargout}] = ft_notification(varargin{:});
elseif isequal(varargin, {'last'})
  % return an answer anyway
  varargout{1} = ft_notification(varargin{:});
else
  ft_notification(varargin{:});
end

