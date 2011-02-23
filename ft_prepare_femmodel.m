function [vol, cfg] = ft_prepare_femmodel(cfg, mri)

% FT_PREPARE_FEMMODEL computes the FEM system matrix.
%
% Use as
%   [vol, cfg] = ft_prepare_femmodel(cfg, seg), or
%   [vol, cfg] = ft_prepare_femmodel(cfg, vol), or
%   [vol, cfg] = ft_prepare_femmodel(cfg)
%
% The configuration can contain
%   cfg.tissue         = [1 2 3 ... N], segmentation value of each tissue type
%   cfg.conductivity   = correspondent conductivity matrix
%   cfg.method         = 'fns', 'simbio'
%

% Copyright (C) 2011, Cristiano Micheli

ft_defaults

if strcmp(cfg.method, 'fns')
  
elseif strcmp(cfg.method, 'simbio')
  error('not yet implemented');
else
  error('unsupported method');
end 