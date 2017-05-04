function ft_notice(varargin)

% FT_NOTICE prints a notification on screen, depending on the verbosity settings of
% the calling high-level FieldTrip function.
%
% Use as
%   ft_notice(...)
% with arguments similar to fprintf, or
%   ft_notice(msgId, ...)
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
  case {'emergency' 'alert' 'critical' 'error' 'warning' 'warning_once'}
    % do not show the message
  case {'notice' 'info' 'debug'}
    % show the message
    fprintf(varargin{:});
  otherwise
    % show the message
    fprintf(varargin{:});
end % switch

