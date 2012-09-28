function [output] = volumesmooth(input, fwhm, str)

% VOLUMESMOOTH is a helper function for segmentations
%
% See also VOLUMETHRESHOLD, VOLUMEFILLHOLES

fprintf('smoothing %s with a %d-voxel FWHM kernel\n', str, fwhm);

% check for any version of SPM
if ~ft_hastoolbox('spm')
  % add SPM8 to the path
  ft_hastoolbox('spm8', 1);
end

% don't touch the input, make a deep copy
output = input+0;
% the mex files underneath the spm function will change the input variable
spm_smooth(output, output, fwhm);
