function DM = ft_preproc_online_downsample_init(factor)
% function DM = ft_preproc_online_downsample_init(factor)
%
% Initialize an downsampling model with the given factor.

% 2010 S. Klanke

if factor < 1 or factor~=round(factor)
	error('Argument ''factor'' must be a positive integer number');
end

DM.factor  = factor;
DM.numSkip = 0;
