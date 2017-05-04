function ft_info(varargin)

% FT_INFO prints an informational message on screen, depending on the verbosity
% settings of the calling high-level FieldTrip function.
%
% Use as
%   ft_info(...)
% with arguments similar to fprintf, or
%   ft_info(msgId, ...)
% with arguments similar to warning.
%
% See also FT_ERROR, FT_WARNING, FT_NOTICE, FT_INFO, FT_DEBUG

% get the configuration from the calling function
try
  cfg = evalin('caller', 'cfg');
catch
  cfg = [];
end

% set the default
cfg.feedback = ft_getopt(cfg, 'feedback', 'notice');

switch cfg.feedback
  case {'emergency' 'alert' 'critical' 'error' 'warning' 'warning_once' 'notice'}
    % do not show the message
  case {'info' 'debug'}
    % show the message
    fprintf(varargin{:});
  otherwise
    % show the message
    fprintf(varargin{:});
end % switch

