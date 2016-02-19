% FT_POSTAMBLE_RANDOMSEED is a helper script that stores the state of the
% random number generator which is used for calling rand/randn/randi functions.
%
% Use as
%   ft_preamble randomseed
%   ... regular code goes here ...
%   ft_postamble randomseed
%
% See also FT_PREAMBLE_RANDOMSEED

if exist('ftFuncRandomseed','var')
  cfg.callinfo.randomseed = ftFuncRandomseed;
end
