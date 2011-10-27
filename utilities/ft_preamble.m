function ft_preamble(cmd, varargin)

% FT_PREAMBLE is a helper function that is included in many of the FieldTrip
% functions and which takes care of some general settings and operations at the
% begin of the function
%
% See also FT_POSTAMBLE

% ideally this would be a script, because the local variables would then be
% shared with the calling function. Instead, this is a function which then
% passes the variables explicitely to another script which is eval'ed.

% the following ensures that these scripts are included as dependencies
% when using the MATLAB compiler
%
%#function ft_preamble_help
%#function ft_preamble_trackconfig
%#function ft_preamble_callinfo
%#function ft_preamble_loadvar

global ft_default

ft_default.preamble = varargin;

if exist(['ft_preamble_' cmd], 'file')
  evalin('caller', ['ft_preamble_' cmd]);
end

ft_default = rmfield(ft_default, 'preamble');
