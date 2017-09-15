function data = chanscale_common(cfg, data)

% CHANSCALE_COMMON applies a scaling to specific channel types
%
% Use as
%   data = chanscale_common(cfg, data)
% where
%   cfg.parameter
%   cfg.magscale
%   cfg.gradscale

cfg.magscale  = ft_getopt(cfg, 'magscale', 1);
cfg.gradscale = ft_getopt(cfg, 'gradscale', 1);

if cfg.magscale~=1 || cfg.gradscale~=1
  error('not yet implemented');
end
