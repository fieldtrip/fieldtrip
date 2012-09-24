function [output] = volumesmooth(input, fwhm, str)

% VOLUMESMOOTH is a helper function for segmentations 

fprintf('smoothing %s with a %d-voxel FWHM kernel\n', str, fwhm);

% don't touch the input, make a deep copy
output = input+0;
% the mex files underneath the spm function will change the input variable
spm_smooth(output, output, fwhm);
