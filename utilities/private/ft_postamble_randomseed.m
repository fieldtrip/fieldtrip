% FT_POSTAMBLE_RANDOMSEED is a helper script that reports the state of the
% random number generator used for calling rand/randn/randi functions.
% See FT_PREAMBLE_RANDOMSEED

if exist('ftFuncRandomseed','var')
  cfg.callinfo.randomseed = ftFuncRandomseed;
end