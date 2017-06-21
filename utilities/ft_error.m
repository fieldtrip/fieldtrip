function varargout = ft_error(varargin)

% FT_ERROR prints an error message on screen, just like the standard ERROR function.
%
% Use as
%   ft_error(...)
% with arguments similar to fprintf, or
%   ft_error(msgId, ...)
% with arguments similar to error.
%
% See also FT_ERROR, FT_WARNING, FT_NOTICE, FT_INFO, FT_DEBUG, ERROR, WARNING

if nargout
  varargout{:} = ft_notification(varargin{:});
elseif isequal(varargin, {'last'})
  % return an answer anyway
  varargout{1} = ft_notification(varargin{:});
else
  ft_notification(varargin{:});
end