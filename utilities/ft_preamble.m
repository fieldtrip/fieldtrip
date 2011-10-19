function ft_preamble(cmd, varargin)

% FT_PREAMBLE is a helper function that is included in many of the FieldTrip
% functions and which takes care of some general settings and operations at the
% begin of the function

% ideally this would be a script, because the local variables would then be
% shared with the calling function. Instead, the variables are explicitely
% copied back and forth between the callers workspace and this one.

global ft_default

ft_default.preamble = varargin;

evalin('caller', ['ft_preamble_' cmd]);

ft_default = rmfield(ft_default, 'preamble');
