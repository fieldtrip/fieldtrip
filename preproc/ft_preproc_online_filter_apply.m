function [FM, xf] = ft_preproc_online_filter_apply(FM, x)
% function [FM, xf] = ft_preproc_online_filter_apply(FM, x)
%
% Passes signal x (channels times samples) through the filter,
% returns updated filter model (delay states!) and filtered signal.
%
% 2010 S. Klanke

[dimX, numX] = size(x);

if dimX > numX
	% explicit algorithm - faster for many channels, 1 sample (=fMRI)
	xf = zeros(size(x));

	for k=1:numX;
		z_old = FM.z;
		z0 = x(:,k) - z_old*FM.A2;
		xf(:,k) = z_old * FM.B2 + z0*FM.B1;
		FM.z(:,2:end) = z_old(:,1:end-1);	
		FM.z(:,1) = z0;
	end
else
	% use built-in Matlab stuff - faster for many samples
	[xf, z] = filter(FM.B,FM.A, x, FM.z',2);
	FM.z = z';
end