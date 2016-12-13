function [data] = ft_transform_ODs(cfg, data)
% FT_TRANSFORM_ODS is outdated, use FT_NIRS_TRANSFORM_ODS instead.

warning('ft_transform_ODs is outdated, please use ft_nirs_transform_ODs instead\n');
[data] = ft_nirs_transform_ODs(cfg, data);

