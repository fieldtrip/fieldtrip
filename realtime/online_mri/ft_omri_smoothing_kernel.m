function [kernX, kernY, kernZ, offsets] = ft_omri_smoothing_kernel(FWHM, voxdims)

% function [kernX, kernY, kernZ, offsets] = ft_omri_smoothing_kernel(FWHM, voxdims)

% Copied from SPM8

if numel(FWHM) == 1
	s = FWHM*ones(1,3);
elseif numel(FWHM) == 3
    s = FHWM;
else
	error 'Invalid size of first parameter';
end

s  = s(:)./voxdims(:);              % voxel anisotropy
s1 = s/sqrt(8*log(2));              % FWHM -> Gaussian parameter

QUA = 3;

x  = round(QUA*s1(1)); x = -x:x; x = spm_smoothkern(s(1),x,1); kernX = single(x/sum(x));
y  = round(QUA*s1(2)); y = -y:y; y = spm_smoothkern(s(2),y,1); kernY = single(y/sum(y));
z  = round(QUA*s1(3)); z = -z:z; z = spm_smoothkern(s(3),z,1); kernZ = single(z/sum(z));

i  = (length(kernX) - 1)/2;
j  = (length(kernY) - 1)/2;
k  = (length(kernZ) - 1)/2;

offsets = [i j k];
