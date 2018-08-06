function [montage, cfg] = ft_prepare_ODtransformation(cfg, data)
% FT_PREPARE_ODTRANSFORMATION is outdated, please use
% FT_NIRS_PREPARE_ODTRANSFORMATION.

warning('ft_prepare_ODtransformation has been renamed to ft_nirs_prepare_ODtransformation\n');
[montage, cfg] = ft_nirs_prepare_ODtransformation(cfg, data);