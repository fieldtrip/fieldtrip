function ft_postamble(cmd, varargin)

% FT_PREAMBLE is a helper function that is included in many of the FieldTrip
% functions and which takes care of some general settings and operations at the
% begin of the function
%
% See also FT_PREAMBLE

% ideally this would be a script, because the local variables would then be
% shared with the calling function. Instead, this is a function which then
% passes the variables explicitely to another script which is eval'ed.

% the following ensures that these scripts are included as dependencies
% when using the MATLAB compiler
%
%#function ft_postamble_trackconfig
%#function ft_postamble_callinfo
%#function ft_postamble_previous
%#function ft_postamble_history
%#function ft_postamble_savevar

global ft_default

% this is a trick to pass the input arguments into the ft_postamble_xxx script
ft_default.postamble = varargin;

if exist(['ft_postamble_' cmd], 'file')
  evalin('caller', ['ft_postamble_' cmd]);
end

if isfield(ft_default, 'postamble')
  % the postamble field should not remain in the ft_default structure
  ft_default = rmfield(ft_default, 'postamble');
end