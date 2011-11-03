% FT_PREAMBLE_TRACKCONFIG is a helper script that calls ft_checkconfig to
% switch the (optional) configuration tracking on. This should be used
% together with FT_POSTAMBLE_TRACKCONFIG.

% otherwise the empty field would end up in the output cfg
global ft_default
ft_default = rmfield(ft_default, 'preamble');

% most fieldtrip functions should allow for configuration tracking, except for
% the functions that take a cfg as input and return a cfg as output
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% the calling ft_preable expects it to be present
ft_default.preamble = {};
